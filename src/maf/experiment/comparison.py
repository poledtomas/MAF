"""Strict compatibility and aligned-bin rebinning from saved cluster statistics."""
import numpy as np
from ..observables import Observable
from ..statistics.errors import jackknife
from ..statistics.aggregation import observable_config
from ..io.metadata import value_at


def definition_key(definition):
    obs = Observable(**definition)
    result = observable_config(obs)
    result.pop("name")
    result.pop("bins")
    # Resolve selection defaults so explicit and implicit inclusive_final agree.
    from dataclasses import asdict
    result["selection"] = asdict(obs._cut)
    cut = result["selection"]
    cut["frame"] = cut["input_frame"] if cut["frame"] == "input" else cut["frame"]
    cut.pop("input_frame")
    cut.pop("cm_rapidity")
    if result["selection"]["pdg"] is not None:
        result["selection"]["pdg"] = sorted(result["selection"]["pdg"])
    return result


def compatibility(result, observable, dataset):
    obs = next(o for o in result.observables if o.name == observable)
    problems = []
    if definition_key(observable_config(obs)) != definition_key(dataset.metadata["definition"]):
        problems.append("Observable, species, acceptance, frame, normalization, units or particle definition differ")
    for key in ("collision_system", "energy", "centrality"):
        values = [value_at(source["metadata"], key) for source in result.sources.values()]
        if not values or any(v is None or v != dataset.metadata[key] for v in values):
            problems.append(f"Missing or incompatible {key}")
    return problems


def rebin(result, observable, edges):
    obs = next(o for o in result.observables if o.name == observable)
    if obs.bins is None or obs.kind in ("mean_pt", "ratio", "mean_mt_minus_m0", "v1", "v2", "vn"):
        raise ValueError("Rebinning requires a count distribution")
    edges = np.asarray(edges, float)
    old = np.asarray(obs.bins)
    indices = []
    for edge in edges:
        match = np.flatnonzero(np.isclose(old,edge,rtol=1e-10,atol=1e-12))
        if len(match) != 1:
            raise ValueError("Cannot split simulation bins or extrapolate: new edges must align")
        indices.append(match[0])
    if len(indices)<2 or np.any(np.diff(indices)<=0):
        raise ValueError("New bin edges must increase")
    offset = sum(o.accumulator_size for o in result.observables[:result.observables.index(obs)])
    a,b = result.arrays()
    a,b = a[:,offset:offset+obs.size], b[:,offset:offset+obs.size]
    if not np.allclose(b,b[:,:1]):
        raise ValueError("Cannot rebin observables with bin-dependent denominators")
    numerators = np.column_stack([a[:,lo:hi].sum(axis=1) for lo,hi in zip(indices[:-1],indices[1:])])
    denominators = np.repeat(b[:,:1],len(indices)-1,axis=1)
    value,error,cov = jackknife(numerators,denominators)
    definition = observable_config(obs)
    if definition.get("scale") is not None and np.asarray(definition["scale"]).size > 1:
        raise ValueError("Rebinning per-bin custom scales is unsupported")
    definition["bins"] = edges.tolist()
    factor = Observable(**definition).normalization()
    return value*factor, error*abs(factor), cov*np.outer(factor,factor)


def compare(result, observable, dataset):
    problems = compatibility(result,observable,dataset)
    if problems:
        raise ValueError("Incompatible comparison: "+"; ".join(problems))
    if np.any(dataset.low[1:] != dataset.high[:-1]):
        raise ValueError("Comparison currently requires contiguous experimental bins")
    edges = np.r_[dataset.low, dataset.high[-1]]
    value,error,covariance = rebin(result,observable,edges)
    with np.errstate(divide="ignore",invalid="ignore"):
        ratio = np.where(dataset.value != 0,value/dataset.value,np.nan)
        model_ratio_error = np.where(dataset.value != 0,error/np.abs(dataset.value),np.nan)
    return {"model": value, "model_stat": error, "model_covariance": covariance,
            "ratio": ratio, "ratio_model_stat": model_ratio_error,
            "data": dataset.value, "data_stat": dataset.stat, "data_systematic": dataset.systematic}


def chi_square(model, data, covariance):
    """User supplies the FULL covariance and its systematics convention."""
    delta = np.asarray(model,float) - np.asarray(data,float)
    cov = np.asarray(covariance,float)
    if delta.ndim != 1 or cov.shape != (len(delta),len(delta)) or not np.isfinite(cov).all() or not np.isfinite(delta).all() or not np.allclose(cov,cov.T):
        raise ValueError("Invalid data or covariance")
    try:
        factor = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Covariance must be positive definite; singular fits need an explicit prescription") from exc
    residual = np.linalg.solve(factor,delta)
    return {"chi2": float(residual@residual), "n_bins": len(delta)}
