"""Rapidity histogram. All bins use [lower, upper)."""
import numpy as np


def calculate(rapidity, bins):
    rapidity = np.asarray(rapidity)
    accepted = np.isfinite(rapidity) & (rapidity >= bins[0]) & (rapidity < bins[-1])
    counts = np.histogram(rapidity[accepted], bins=bins)[0]
    return counts, 1.0


def normalization(bins):
    return 1.0 / np.diff(bins)
