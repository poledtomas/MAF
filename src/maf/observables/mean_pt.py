"""Particle-weighted mean pT: ensemble sum(pT) / ensemble particle count."""
import numpy as np


def calculate(pt):
    pt = np.asarray(pt)
    return [pt.sum()], len(pt)
