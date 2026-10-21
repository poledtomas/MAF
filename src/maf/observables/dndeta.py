"""Pseudorapidity histogram; the main analysis selects charged particles."""
import numpy as np


def calculate(eta, bins):
    eta = np.asarray(eta)
    accepted = np.isfinite(eta) & (eta >= bins[0]) & (eta < bins[-1])
    counts = np.histogram(eta[accepted], bins=bins)[0]
    return counts, 1.0


def normalization(bins):
    return 1.0 / np.diff(bins)
