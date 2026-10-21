"""Independent y, eta and pT acceptance; lower inclusive, upper exclusive."""
from dataclasses import dataclass
import numpy as np


def all_charged_alias(pdg):
    """Accept pdg: all_charged and pdg: [all_charged]."""
    return ((isinstance(pdg,str) and pdg=='all_charged') or
            (isinstance(pdg,(list,tuple)) and len(pdg)==1 and
             isinstance(pdg[0],str) and pdg[0]=='all_charged'))


def normalize_selection(settings):
    """Store the alias canonically as charged: true, with no species filter."""
    result=dict(settings)
    pdg=result.get('pdg')
    if all_charged_alias(pdg):
        if result.get('charged') is not None and not isinstance(result['charged'],(bool,np.bool_)):
            raise ValueError('charged must be a boolean')
        if result.get('charged') is not None and not result['charged']:
            raise ValueError('all_charged conflicts with charged: false')
        result.pop('pdg')
        result['charged']=True
    elif isinstance(pdg,(list,tuple)) and any(isinstance(p,str) and p=='all_charged' for p in pdg):
        raise ValueError('all_charged must be used alone, without numeric PDG IDs')
    return result


@dataclass
class Selection:
    pdg: object = None
    charged: object = None
    y: object = None
    eta: object = None
    pt: object = None
    frame: str = "input"
    input_frame: str = "unknown"
    cm_rapidity: object = None
    particle_definition: str = "inclusive_final"
    boundary: str = "half_open"

    def __post_init__(self):
        normalized=normalize_selection({'pdg':self.pdg,'charged':self.charged})
        self.pdg=normalized.get('pdg')
        self.charged=normalized.get('charged')
        if self.boundary not in ("half_open", "open"):
            raise ValueError("boundary must be half_open or open")
        if self.charged is not None and not isinstance(self.charged, (bool, np.bool_)):
            raise ValueError("charged must be a boolean")
        if self.pdg is not None:
            pdg = np.asarray(self.pdg)
            if pdg.ndim != 1 or pdg.dtype.kind not in "iu" or pdg.dtype.kind == "b":
                raise ValueError("pdg must be a list of integer PDG IDs")
        if self.particle_definition != "inclusive_final":
            raise ValueError("Only inclusive_final particles are supported; no inferred feed-down corrections")
        for key in ("y", "eta", "pt"):
            value = getattr(self, key)
            if value is not None and (len(value) != 2 or not np.isfinite(value).all() or value[0] >= value[1]):
                raise ValueError(f"Invalid {key} acceptance")
        if self.pt is not None and self.pt[0] < 0:
            raise ValueError("pT acceptance must be nonnegative")
        if self.frame not in ("input", "lab", "cm") or self.input_frame not in ("unknown", "lab", "cm"):
            raise ValueError("Frame must be input, lab or cm")
        if self.frame != "input" and self.input_frame == "unknown":
            raise ValueError("An explicit input_frame is required for frame selection")
        if self.frame not in ("input", self.input_frame) and (self.cm_rapidity is None or not np.isfinite(self.cm_rapidity)):
            raise ValueError("Frame conversion requires finite cm_rapidity in the lab frame")
        if self.eta is not None and self.frame not in ("input", self.input_frame):
            raise ValueError("Eta cuts across frames require a full momentum boost; use the input frame")

    def rapidity(self, p):
        y = p.rapidity
        if self.frame == "cm" and self.input_frame == "lab":
            return y - self.cm_rapidity
        if self.frame == "lab" and self.input_frame == "cm":
            return y + self.cm_rapidity
        return y

    def mask(self, p):
        mask = np.ones(len(p), dtype=bool)
        if self.pdg is not None:
            mask &= np.isin(p.pdg, self.pdg)
        if self.charged is not None:
            mask &= (p.charge != 0) if self.charged else (p.charge == 0)
        for key, values in (("y", self.rapidity(p)), ("eta", p.eta), ("pt", p.pt)):
            bounds = getattr(self, key)
            if bounds is not None:
                lower = values >= bounds[0] if self.boundary == "half_open" else values > bounds[0]
                mask &= np.isfinite(values) & lower & (values < bounds[1])
        return mask
