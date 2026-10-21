"""Equal-subevent event-plane resolution, following Poskanzer and Voloshin."""
import numpy as np


def resolution(chi):
    try:
        from scipy.special import i0e, i1e
    except ImportError as exc:
        raise ImportError("Event-plane flow requires maf[flow] (SciPy)") from exc
    chi = np.asarray(chi, float)
    argument = chi**2 / 4
    return np.sqrt(np.pi) / (2*np.sqrt(2)) * chi * (i0e(argument) + i1e(argument))


def full_resolution(correlation):
    """Invert R(chi_sub)=sqrt(<cos n(Psi_A-Psi_B)>), then use sqrt(2)*chi_sub."""
    correlation = np.asarray(correlation, float)
    valid = np.isfinite(correlation) & (correlation > 0) & (correlation <= 1 + 1e-12)
    target = np.sqrt(np.clip(correlation, 0, 1))
    low, high = np.zeros_like(target), np.full_like(target, 1e6)
    for _ in range(65):
        middle = (low+high)/2
        below = resolution(middle) < target
        low, high = np.where(below,middle,low), np.where(below,high,middle)
    result = resolution(np.sqrt(2)*(low+high)/2)
    result = np.where(correlation >= 1, 1.0, result)
    return np.where(valid, result, np.nan)
