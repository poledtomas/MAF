"""Known-plane v_n for any positive integer harmonic, versus y, eta or pT."""
import numpy as np


def calculate(px, py, coordinate, bins, plane_angle, harmonic, weighting="particle"):
    if type(harmonic) is not int or harmonic < 1:
        raise ValueError("Flow harmonic must be a positive integer")
    if not np.isscalar(plane_angle) or np.asarray(plane_angle).dtype.kind not in "fiu" or not np.isfinite(plane_angle):
        raise ValueError("Reaction-plane angle must be one finite number in radians")
    px, py, coordinate = np.broadcast_arrays(px, py, coordinate)
    accepted = np.isfinite(coordinate) & np.isfinite(px) & np.isfinite(py) & (np.hypot(px, py) > 0)
    accepted &= (coordinate >= bins[0]) & (coordinate < bins[-1])
    phi = np.arctan2(py[accepted], px[accepted])
    numerator = np.histogram(coordinate[accepted], bins=bins, weights=np.cos(harmonic*(phi-plane_angle)))[0]
    denominator = np.histogram(coordinate[accepted], bins=bins)[0]
    if weighting == "event":
        numerator = np.divide(numerator, denominator, out=np.zeros(len(bins)-1), where=denominator > 0)
        denominator = denominator > 0
    elif weighting != "particle":
        raise ValueError("Flow weighting must be event or particle")
    return numerator, denominator
