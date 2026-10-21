"""Uncertainties from independent clusters, preserving paired statistics."""
import numpy as np


def ratio_estimate(numerator, denominator):
    a, b = np.asarray(numerator, float), np.asarray(denominator, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(a, b, out=np.full(np.broadcast_shapes(a.shape, b.shape), np.nan), where=b != 0)


def jackknife(numerator, denominator):
    """Paired delete-one-unit uncertainty for a ratio of sums."""
    return jackknife_derived(numerator, denominator, ratio_estimate)


def jackknife_derived(numerator, denominator, estimator):
    """Estimator takes aggregate numerator/denominator arrays with features last."""
    a, b = np.asarray(numerator,float), np.asarray(denominator,float)
    if a.ndim == 1:
        a, b = a[:,None], b[:,None]
    total_a, total_b = a.sum(axis=0), b.sum(axis=0)
    estimate = np.atleast_1d(estimator(total_a,total_b))
    n, k = len(a), len(estimate)
    covariance = np.full((k,k),np.nan)
    if n < 2:
        return estimate, np.full(k,np.nan), covariance
    replicas = estimator(total_a-a,total_b-b)
    valid = np.isfinite(replicas).all(axis=0) & np.isfinite(estimate)
    centered = replicas[:,valid] - replicas[:,valid].mean(axis=0)
    covariance[np.ix_(valid,valid)] = (n-1)/n * centered.T@centered
    return estimate, np.sqrt(np.maximum(0,np.diag(covariance))), covariance


def bootstrap(numerator, denominator, *, samples=1000, seed=0):
    """Paired independent-cluster bootstrap; returns replicas, never particle resampling."""
    a, b = np.asarray(numerator, float), np.asarray(denominator, float)
    if len(a) < 2 or samples < 2:
        raise ValueError("Bootstrap requires at least two independent units and replicas")
    rng = np.random.default_rng(seed)
    result = []
    for _ in range(samples):
        indices = rng.integers(len(a), size=len(a))
        result.append(ratio_estimate(a[indices].sum(axis=0), b[indices].sum(axis=0)))
    return np.asarray(result)
