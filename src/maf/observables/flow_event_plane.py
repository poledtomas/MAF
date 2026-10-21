"""pT-weighted EP, leave-one-particle-out autocorrelation removal and parity subevents."""
import numpy as np
from ..statistics.errors import ratio_estimate
from .flow_resolution import full_resolution


def calculate(phi, pt, coordinate, reference, particles, bins, harmonic, weighting="event", empty_bins="exclude"):
    if type(harmonic) is not int or harmonic < 1:
        raise ValueError("Flow harmonic must be a positive integer")
    if weighting not in ("event", "particle") or empty_bins not in ("exclude", "zero"):
        raise ValueError("Invalid event-plane weighting or empty-bin policy")
    if empty_bins == "zero" and weighting != "event":
        raise ValueError("empty_bins: zero requires event weighting")
    numerator, denominator = np.zeros(len(bins)), np.zeros(len(bins))
    vector = pt * np.exp(1j*harmonic*phi)
    indices = np.arange(len(phi))
    a = vector[reference & (indices % 2 == 0)].sum()
    b = vector[reference & (indices % 2 == 1)].sum()
    tolerance = 1e-12 * pt[reference].sum()
    if abs(a) <= tolerance or abs(b) <= tolerance:
        return numerator, denominator
    # Resolution uses all eligible events, independently of occupancy of the first pT bin.
    numerator[-1] = np.real(a*np.conj(b)) / (abs(a)*abs(b))
    denominator[-1] = 1.0
    plane = vector[reference].sum() - np.where(reference, vector, 0)
    accepted = particles & np.isfinite(coordinate) & (abs(plane) > tolerance)
    accepted &= (coordinate >= bins[0]) & (coordinate < bins[-1])
    cosine = np.real(np.exp(1j*harmonic*phi[accepted]) * np.conj(plane[accepted])/abs(plane[accepted]))
    numerator[:-1] = np.histogram(coordinate[accepted],bins=bins,weights=cosine)[0]
    denominator[:-1] = np.histogram(coordinate[accepted],bins=bins)[0]
    if weighting == "event":
        numerator[:-1] = np.divide(numerator[:-1],denominator[:-1],out=np.zeros(len(bins)-1),where=denominator[:-1]>0)
        denominator[:-1] = 1.0 if empty_bins == "zero" else denominator[:-1] > 0
    return numerator, denominator


def finalize(numerator, denominator):
    ratios = ratio_estimate(numerator, denominator)
    correction = full_resolution(ratios[..., -1])
    with np.errstate(divide="ignore",invalid="ignore"):
        return ratios[..., :-1] / correction[..., None]
