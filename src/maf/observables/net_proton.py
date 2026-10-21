"""Net protons counted within each event to preserve p/pbar covariance."""
import numpy as np


def calculate(rapidity, pdg, bins):
    rapidity, pdg = np.asarray(rapidity), np.asarray(pdg)
    accepted = (np.abs(pdg) == 2212) & np.isfinite(rapidity)
    accepted &= (rapidity >= bins[0]) & (rapidity < bins[-1])
    signs = np.where(pdg[accepted] == 2212, 1.0, -1.0)
    counts = np.histogram(rapidity[accepted], bins=bins, weights=signs)[0]
    return counts, 1.0


def normalization(bins):
    return 1.0 / np.diff(bins)
