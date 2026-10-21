"""Transverse momentum spectra. Momenta and bin edges are in GeV."""
import numpy as np


def calculate(pt, bins):
    pt = np.asarray(pt)
    accepted = np.isfinite(pt) & (pt >= bins[0]) & (pt < bins[-1])
    counts = np.histogram(pt[accepted], bins=bins)[0]
    return counts, 1.0


def normalization(bins, invariant=False):
    edges = np.asarray(bins)
    if invariant:
        # Exact finite-bin annular area. Rapidity width is applied by the driver.
        return 1.0 / (np.pi * (edges[1:]**2 - edges[:-1]**2))
    return 1.0 / np.diff(edges)
