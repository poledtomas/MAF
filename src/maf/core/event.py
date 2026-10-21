"""Input-format-independent particle and event containers; GeV, c=1."""
from dataclasses import dataclass, field
from typing import Optional
import numpy as np


def longitudinal_rapidity(energy, pz):
    energy, pz = np.broadcast_arrays(energy, pz)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 0.5 * (np.log(energy + pz) - np.log(energy - pz))


@dataclass
class Particles:
    energy: np.ndarray
    px: np.ndarray
    py: np.ndarray
    pz: np.ndarray
    mass: np.ndarray
    pdg: np.ndarray
    charge: np.ndarray
    ids: Optional[np.ndarray] = None
    extras: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in ("energy", "px", "py", "pz", "mass", "pdg", "charge"):
            raw = np.asarray(getattr(self, name))
            if name in ("pdg", "charge") and (not np.isfinite(raw).all() or np.any(raw != np.floor(raw))):
                raise ValueError(f"{name} must contain integers")
            setattr(self, name, np.asarray(raw, dtype=np.int64 if name in ("pdg", "charge") else float))
        arrays = [getattr(self, n) for n in ("energy", "px", "py", "pz", "mass", "pdg", "charge")]
        if any(a.ndim != 1 or len(a) != len(self.energy) for a in arrays):
            raise ValueError("Particle arrays must be one-dimensional with equal lengths")
        if any(not np.isfinite(a).all() for a in arrays):
            raise ValueError("Nonfinite particle data")
        if np.any(self.energy < 0) or np.any(self.mass < 0):
            raise ValueError("Negative energy or mass")
        if self.ids is not None:
            self.ids = np.asarray(self.ids, dtype=np.int64)
            if self.ids.shape != self.energy.shape:
                raise ValueError("Particle IDs must match particle count")

    def __len__(self):
        return len(self.energy)

    @property
    def pt(self):
        return np.hypot(self.px, self.py)

    @property
    def rapidity(self):
        return longitudinal_rapidity(self.energy, self.pz)

    @property
    def eta(self):
        # Stable at large |pz|; a particle at rest has undefined eta.
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.arcsinh(self.pz / self.pt)

    def mt(self, mass=None):
        return np.hypot(self.mass if mass is None else mass, self.pt)

    def validate_shell(self, rtol=5e-4, atol=1e-6):
        expected = np.sqrt(self.mass**2 + self.px**2 + self.py**2 + self.pz**2)
        if not np.allclose(self.energy, expected, rtol=rtol, atol=atol):
            raise ValueError("Energy and momentum violate the mass shell within configured tolerance")


@dataclass
class Event:
    particles: Particles
    event_id: int
    ensemble_id: int = 0
    weight: float = 1.0
    hydro_id: Optional[str] = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if not np.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("Event weights must be finite and positive")
