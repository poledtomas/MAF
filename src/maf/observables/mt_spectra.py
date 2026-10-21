"""Transverse mass spectra with an explicit reference mass in GeV."""
import numpy as np


def calculate(pt, mass, bins, subtract_mass=False):
    mt = np.hypot(mass, np.asarray(pt))
    values = mt - mass if subtract_mass else mt
    accepted = np.isfinite(values) & (values >= bins[0]) & (values < bins[-1])
    counts = np.histogram(values[accepted], bins=bins)[0]
    return counts, 1.0


def normalization(bins, invariant=False):
    edges = np.asarray(bins)
    if invariant:
        return 1.0 / (np.pi * (edges[1:]**2 - edges[:-1]**2))
    return 1.0 / np.diff(edges)
