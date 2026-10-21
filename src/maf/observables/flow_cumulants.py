"""Unit-weight direct Q cumulants, including arbitrary POI/reference overlap.

Bilandzic et al., arXiv:1010.0233v2, Eqs. (16), (18), (24), (25),
(28), (32), (34), (35).
"""
import numpy as np
from ..statistics.errors import ratio_estimate


def calculate(phi, coordinate, reference, particles, bins, harmonic, four=False):
    nbin = len(bins)-1
    numerator = np.zeros(2+2*nbin if four else 1+nbin)
    denominator = np.zeros_like(numerator)
    multiplicity = int(np.count_nonzero(reference))
    if multiplicity == 0:
        return numerator, denominator
    phase = np.exp(1j*harmonic*phi)
    q = phase[reference].sum()
    q2 = (phase[reference]**2).sum()
    q_squared = abs(q)**2
    if multiplicity >= 2:
        numerator[0] = q_squared-multiplicity
        denominator[0] = multiplicity*(multiplicity-1)
    offset = 2 if four else 1
    if four and multiplicity >= 4:
        numerator[1] = q_squared**2 + abs(q2)**2 - 2*np.real(q2*np.conj(q)**2) - 4*(multiplicity-2)*q_squared + 2*multiplicity*(multiplicity-3)
        denominator[1] = multiplicity*(multiplicity-1)*(multiplicity-2)*(multiplicity-3)
    for index, (low, high) in enumerate(zip(bins[:-1],bins[1:])):
        selected = particles & np.isfinite(coordinate) & (coordinate >= low) & (coordinate < high)
        overlap = selected & reference
        count = int(np.count_nonzero(selected))
        overlap_count = int(np.count_nonzero(overlap))
        p = phase[selected].sum()
        overlap_q = phase[overlap].sum()
        overlap_q2 = (phase[overlap]**2).sum()
        # Eqs. (24), (28): exclude only particles actually shared with RFPs.
        pair_weight = count*multiplicity-overlap_count
        if pair_weight > 0:
            numerator[offset+index] = np.real(p*np.conj(q))-overlap_count
            denominator[offset+index] = pair_weight
        if four and multiplicity >= 3:
            four_weight = (count*multiplicity-3*overlap_count)*(multiplicity-1)*(multiplicity-2)
            if four_weight > 0:
                # Eq. (32). p_n includes all POIs; q_n and q_2n here include
                # only the POI/RFP intersection. They need not equal p_n.
                value = p*q*np.conj(q)**2 - overlap_q2*np.conj(q)**2 - p*q*np.conj(q2)
                value -= 2*(multiplicity-1)*p*np.conj(q) + 2*overlap_count*q_squared
                value += 7*overlap_q*np.conj(q) - q*np.conj(overlap_q) + overlap_q2*np.conj(q2)
                numerator[offset+nbin+index] = np.real(value) + 2*overlap_count*(multiplicity-3)
                denominator[offset+nbin+index] = four_weight
    return numerator, denominator


def finalize(numerator, denominator, four=False):
    """Return differential v_n{2/4}, not the raw four-particle correlation.

    Inputs are ensemble tuple sums. Eqs. (12), (34), (35) require products
    of ensemble averages, not averages of eventwise cumulants or v_n.
    """
    ratios = ratio_estimate(numerator, denominator)
    c2 = ratios[..., 0]
    with np.errstate(divide="ignore",invalid="ignore"):
        if not four:
            return np.where(c2[...,None]>0,ratios[...,1:]/np.sqrt(c2[...,None]),np.nan)
        nbin = (ratios.shape[-1]-2)//2
        c4 = ratios[...,1] - 2*c2**2
        d4 = ratios[...,2+nbin:] - 2*ratios[...,2:2+nbin]*c2[...,None]
        return np.where(c4[...,None]<0,-d4/(-c4[...,None])**0.75,np.nan)
