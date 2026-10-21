"""Shared adapter: validate settings, select particles and call plain functions."""
from dataclasses import dataclass, field, replace
from typing import ClassVar
import numpy as np
from ..core.selection import Selection, normalize_selection
from . import dndy, dndeta, pt_spectra, mt_spectra, yields, mean_pt, mean_mt, ratios, net_proton, flow_settings

_BUILTINS = {"yield", "charged_multiplicity", "mean_pt", "mean_mt_minus_m0", "ratio",
             "dndy", "dndeta", "net_proton", "pt", "invariant_pt", "mt", "mt_minus_m0", "invariant_mt", "v1", "v2", "vn"}

_REGISTRY = {}


def register_observable(name, evaluator):
    """Register evaluator(observable, event) -> (numerator, denominator).

    Vectors must match Observable.size. Declare normalization through scale.
    """
    if name in _REGISTRY or name in _BUILTINS:
        raise ValueError(f"Observable already registered: {name}")
    _REGISTRY[name] = evaluator


@dataclass
class Observable:
    name: str
    kind: str
    selection: dict = field(default_factory=dict)
    bins: object = None
    mass: object = None
    denominator_pdg: object = None
    density: object = None
    units: str = ""
    scale: object = None
    axis: str = "y"
    flow: object = None
    harmonic: object = None
    _cut: Selection = field(init=False, repr=False)
    _reference: object = field(init=False, default=None, repr=False)
    required_fields: ClassVar[tuple] = ("energy", "px", "py", "pz", "mass", "pdg", "charge")

    def __post_init__(self):
        if self.kind not in _REGISTRY and self.kind not in _BUILTINS:
            raise ValueError(f"Unknown observable: {self.kind}")
        if not self.name or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in self.name):
            raise ValueError("Observable name must be a safe filename component")
        self.selection = normalize_selection(self.selection)
        if self.kind in ("dndeta", "charged_multiplicity"):
            if self.selection.get("charged") is False:
                raise ValueError("Charged observables cannot select neutral particles")
            self.selection["charged"] = True
        self._cut = Selection(**self.selection)
        if self.kind == "dndeta" and self._cut.frame not in ("input", self._cut.input_frame):
            raise ValueError("dN/deta requires the input momentum frame")
        if self.bins is not None and self.kind in ("yield", "charged_multiplicity", "mean_pt", "mean_mt_minus_m0", "ratio"):
            raise ValueError("Scalar observables do not accept histogram bins")
        if self.bins is not None:
            self.bins = np.asarray(self.bins, dtype=float).tolist()
            edges = np.asarray(self.bins)
            if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
                raise ValueError("Bin edges must be finite and strictly increasing")
        if self.kind in ("dndy", "net_proton", "dndeta", "pt", "invariant_pt", "mt", "mt_minus_m0", "invariant_mt", "v1", "v2", "vn") and self.bins is None:
            raise ValueError("Distribution requires bin edges")
        if self.kind in ("v1", "v2", "vn"):
            flow_settings.validate(self)
        elif self.flow is not None or self.axis != "y" or self.harmonic is not None:
            raise ValueError("axis, harmonic and flow settings apply only to v1/v2/vn")
        if self.kind in ("pt", "invariant_pt", "mt", "mt_minus_m0", "invariant_mt") and self.bins[0] < 0:
            raise ValueError("Spectrum bins must be nonnegative")
        if self.kind in ("invariant_pt", "invariant_mt") and self._cut.y is None:
            raise ValueError("Invariant spectra require an explicit rapidity window")
        if self.kind in ("mt", "mt_minus_m0", "mean_mt_minus_m0", "invariant_mt") and (self.mass is None or not np.isfinite(self.mass) or self.mass <= 0):
            raise ValueError("mT observables require an explicit positive rest mass in GeV")
        if self.mass is not None and self.kind in ("mt", "mt_minus_m0", "mean_mt_minus_m0", "invariant_mt"):
            if self._cut.pdg is None or len(set(abs(p) for p in self._cut.pdg)) != 1:
                raise ValueError("mT with a reference mass requires a single species (optionally including antiparticles)")
        if self.kind == "ratio" and (not self.denominator_pdg or self._cut.pdg is None):
            raise ValueError("Ratio requires numerator selection.pdg and denominator_pdg")
        if self.kind == "net_proton" and self._cut.pdg is not None:
            raise ValueError("Net protons select both p and pbar internally; omit selection.pdg")
        if self.density is not None and self.kind in ("dndy", "dndeta", "net_proton", "invariant_pt", "invariant_mt", "mean_pt", "mean_mt_minus_m0", "ratio", "v1", "v2", "vn"):
            raise ValueError("Additional density normalization is not defined for this observable")
        if self.kind in ("mt", "invariant_mt") and self.bins[0] < self.mass:
            raise ValueError("mT bin edges cannot extend below the reference rest mass")
        if self.density is not None and self.density not in ("y", "eta"):
            raise ValueError("density must be y or eta")
        if self.density is not None and getattr(self._cut, self.density) is None:
            raise ValueError("Density normalization requires the corresponding acceptance window")
        # Truncated bins cannot be normalized as if their full width were used.
        axis = {"dndy": "y", "net_proton": "y", "dndeta": "eta", "pt": "pt", "invariant_pt": "pt"}.get(self.kind)
        if self.kind in ("v1", "v2", "vn"):
            axis = self.axis
        if axis and self.bins is not None and getattr(self._cut, axis) is not None:
            bounds = getattr(self._cut, axis)
            if self.bins[0] < bounds[0] or self.bins[-1] > bounds[1]:
                raise ValueError("Histogram bins extend outside the acceptance on that axis")
        if self.kind in ("mt", "invariant_mt", "mt_minus_m0") and self._cut.pt is not None:
            transformed = np.hypot(self.mass, np.asarray(self._cut.pt))
            if self.kind == "mt_minus_m0":
                transformed -= self.mass
            if self.bins[0] < transformed[0] or self.bins[-1] > transformed[1]:
                raise ValueError("mT bins extend outside the transformed pT acceptance")
        if not self.units:
            self.units = self.default_units()
        if self.scale is not None and (np.asarray(self.scale).size not in (1, self.size) or not np.isfinite(self.scale).all()):
            raise ValueError("Invalid observable scale")

    @property
    def size(self):
        return len(self.bins)-1 if self.bins is not None else 1

    @property
    def accumulator_labels(self):
        if self.kind in ("v1","v2","vn"):
            return flow_settings.labels(self)
        return [f"bin_{i}" for i in range(self.size)]

    @property
    def accumulator_size(self):
        return len(self.accumulator_labels)

    def estimate(self,numerator,denominator):
        if self.kind in ("v1","v2","vn"):
            return flow_settings.finalize(self,numerator,denominator)
        from ..statistics.errors import ratio_estimate
        return ratio_estimate(numerator,denominator)

    def diagnostics(self,numerator,denominator):
        if self.kind in ("v1","v2","vn"):
            return flow_settings.diagnostics(self,numerator,denominator)
        return {}

    def default_units(self):
        if self.kind in ("mean_pt", "mean_mt_minus_m0"):
            return "GeV"
        if self.kind in ("invariant_pt", "invariant_mt"):
            return "GeV^-2"
        if self.kind in ("pt", "mt", "mt_minus_m0"):
            return "GeV^-1"
        return "1"

    def normalization(self):
        if self.kind in ("v1", "v2", "vn"):
            factor = np.ones(self.size)  # flow is a mean cosine, not a density
        elif self.kind in ("pt", "invariant_pt"):
            factor = pt_spectra.normalization(self.bins, invariant=self.kind == "invariant_pt")
        elif self.kind in ("mt", "mt_minus_m0", "invariant_mt"):
            factor = mt_spectra.normalization(self.bins, invariant=self.kind == "invariant_mt")
        elif self.kind == "dndy":
            factor = dndy.normalization(self.bins)
        elif self.kind == "dndeta":
            factor = dndeta.normalization(self.bins)
        elif self.kind == "net_proton":
            factor = net_proton.normalization(self.bins)
        elif self.bins is not None:
            factor = 1.0 / np.diff(self.bins)
        else:
            factor = np.ones(self.size)
        density = "y" if self.kind in ("invariant_pt", "invariant_mt") else self.density
        if density:
            bounds = getattr(self._cut, density)
            factor /= bounds[1] - bounds[0]
        if self.scale is not None:
            factor *= self.scale
        return factor

    def evaluate(self, event):
        if self.kind in _REGISTRY:
            numerator, denominator = _REGISTRY[self.kind](self, event)
        else:
            numerator, denominator = calculate_event(self, event)
        numerator = np.atleast_1d(np.asarray(numerator, dtype=float))
        denominator = np.broadcast_to(np.asarray(denominator, dtype=float), (self.accumulator_size,)).copy()
        if numerator.shape != (self.accumulator_size,) or not np.isfinite(numerator).all() or not np.isfinite(denominator).all() or np.any(denominator < 0):
            raise ValueError(f"Invalid event statistics from {self.name}")
        return numerator, denominator


def calculate_event(observable, event):
    """Apply common cuts once, then call the selected small calculation module."""
    obs = observable
    p = event.particles
    mask = obs._cut.mask(p)
    if obs.kind in ("yield", "charged_multiplicity"):
        return yields.calculate(np.count_nonzero(mask))
    if obs.kind == "mean_pt":
        return mean_pt.calculate(p.pt[mask])
    if obs.kind == "mean_mt_minus_m0":
        return mean_mt.calculate(p.pt[mask], obs.mass)
    if obs.kind == "ratio":
        denominator_cut = replace(obs._cut, pdg=obs.denominator_pdg)
        return ratios.calculate(np.count_nonzero(mask), np.count_nonzero(denominator_cut.mask(p)))
    if obs.kind == "dndy":
        return dndy.calculate(obs._cut.rapidity(p)[mask], obs.bins)
    if obs.kind == "dndeta":
        return dndeta.calculate(p.eta[mask], obs.bins)
    if obs.kind == "net_proton":
        return net_proton.calculate(obs._cut.rapidity(p)[mask], p.pdg[mask], obs.bins)
    if obs.kind in ("pt", "invariant_pt"):
        return pt_spectra.calculate(p.pt[mask], obs.bins)
    if obs.kind in ("mt", "mt_minus_m0", "invariant_mt"):
        return mt_spectra.calculate(p.pt[mask], obs.mass, obs.bins, subtract_mass=obs.kind == "mt_minus_m0")
    if obs.kind in ("v1", "v2", "vn"):
        return flow_settings.calculate(obs,event,mask)
    raise ValueError(f"Unknown observable: {obs.kind}")
