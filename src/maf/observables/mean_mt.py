"""Particle-weighted mean mT minus the reference rest mass."""
import numpy as np


def calculate(pt, mass):
    mt_minus_mass = np.hypot(mass, np.asarray(pt)) - mass
    return [mt_minus_mass.sum()], len(mt_minus_mass)
