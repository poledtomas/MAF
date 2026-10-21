"""Shared flow configuration and array selection; calculation modules stay independent."""
import numpy as np
from ..core.selection import Selection, normalize_selection
from ..statistics.errors import ratio_estimate
from . import flow_reaction_plane, flow_event_plane, flow_cumulants


def validate(obs):
    if obs.kind in ("v1", "v2"):
        fixed = int(obs.kind[1:])
        if obs.harmonic is None:
            obs.harmonic = fixed
        elif type(obs.harmonic) is not int or obs.harmonic != fixed:
            raise ValueError("harmonic must agree with the v1/v2 kind; use kind: vn for other orders")
    if type(obs.harmonic) is not int or obs.harmonic < 1:
        raise ValueError("Flow harmonic must be a positive integer (kind: vn requires harmonic)")
    if obs.axis not in ("y", "eta", "pt"):
        raise ValueError("Flow axis must be y, eta or pt")
    if not isinstance(obs.flow,dict):
        raise ValueError("v1/v2/vn require an explicit flow.method and reference definition")
    obs.flow = dict(obs.flow)
    method = obs.flow.get("method")
    if method not in ("reaction_plane","event_plane","cumulant_2","cumulant_4"):
        raise ValueError("flow.method must be reaction_plane, event_plane, cumulant_2 or cumulant_4")
    allowed = {"method","weighting"} | ({"angle","metadata_key"} if method=="reaction_plane" else {"reference"})
    if method=="event_plane":
        allowed.update({"subevents", "empty_bins"})
    if set(obs.flow)-allowed:
        raise ValueError("Unknown settings for the chosen flow method")
    if method=="reaction_plane":
        if ("angle" in obs.flow)==("metadata_key" in obs.flow):
            raise ValueError("Specify exactly one flow.angle or flow.metadata_key")
        if "angle" in obs.flow:
            angle=obs.flow["angle"]
            if not np.isscalar(angle) or np.asarray(angle).dtype.kind not in "fiu" or not np.isfinite(angle):
                raise ValueError("flow.angle must be finite, in radians")
        if "metadata_key" in obs.flow and (not isinstance(obs.flow["metadata_key"],str) or not obs.flow["metadata_key"]):
            raise ValueError("flow.metadata_key must be a nonempty event metadata key")
        obs.flow.setdefault("weighting","particle")
    else:
        if not isinstance(obs.flow.get("reference"),dict):
            raise ValueError("flow.reference must explicitly define reference particle acceptance")
        obs.flow["reference"] = normalize_selection(obs.flow["reference"])
        obs._reference = Selection(**obs.flow["reference"])
        obs.flow.setdefault("weighting","event" if method=="event_plane" else "tuple")
    expected=("tuple",) if method.startswith("cumulant") else ("event","particle")
    if obs.flow["weighting"] not in expected:
        raise ValueError("Invalid weighting for this flow method")
    if method=="event_plane":
        obs.flow.setdefault("subevents","alternating")
        if obs.flow["subevents"]!="alternating":
            raise ValueError("Only explicitly alternating equal-subevent EP is implemented")
        obs.flow.setdefault("empty_bins", "exclude")
        if obs.flow["empty_bins"] not in ("exclude", "zero"):
            raise ValueError("flow.empty_bins must be exclude or zero")
        if obs.flow["empty_bins"] == "zero" and obs.flow["weighting"] != "event":
            raise ValueError("flow.empty_bins: zero requires event weighting")
        from .flow_resolution import resolution
        resolution(0.0)  # fail at configuration time when optional SciPy is absent
    if obs.axis=="pt" and obs.bins[0]<0:
        raise ValueError("Flow pT bins must be nonnegative")
    if obs.axis=="eta" and obs._cut.frame not in ("input",obs._cut.input_frame):
        raise ValueError("Flow versus eta requires the input momentum frame")


def labels(obs):
    bins=[f"bin_{i}" for i in range(obs.size)]
    method=obs.flow["method"]
    if method=="event_plane":
        return bins+["subevent_correlation"]
    if method=="cumulant_2":
        return ["reference_2"]+["differential_2_"+b for b in bins]
    if method=="cumulant_4":
        return ["reference_2","reference_4"]+["differential_2_"+b for b in bins]+["differential_4_"+b for b in bins]
    return bins


def calculate(obs,event,selected):
    p=event.particles
    coordinate={"y":obs._cut.rapidity(p),"eta":p.eta,"pt":p.pt}[obs.axis]
    method=obs.flow["method"]
    selected=selected & (p.pt>0)
    if method=="reaction_plane":
        plane=obs.flow.get("angle")
        if "metadata_key" in obs.flow:
            key=obs.flow["metadata_key"]
            if key not in event.metadata:
                raise ValueError(f"Missing reaction-plane angle in event metadata: {key}")
            plane=event.metadata[key]
        return flow_reaction_plane.calculate(p.px[selected],p.py[selected],coordinate[selected],
                                             obs.bins,plane,obs.harmonic,obs.flow["weighting"])
    phi=np.arctan2(p.py,p.px)
    reference=obs._reference.mask(p) & (p.pt>0)
    harmonic=obs.harmonic
    if method=="event_plane":
        return flow_event_plane.calculate(phi,p.pt,coordinate,reference,selected,obs.bins,harmonic,obs.flow["weighting"],obs.flow["empty_bins"])
    selected &= np.isfinite(coordinate) & (coordinate>=obs.bins[0]) & (coordinate<obs.bins[-1])
    return flow_cumulants.calculate(phi,coordinate,reference,selected,obs.bins,harmonic,four=method=="cumulant_4")


def finalize(obs,numerator,denominator):
    method=obs.flow["method"]
    if method=="event_plane":
        return flow_event_plane.finalize(numerator,denominator)
    if method.startswith("cumulant"):
        return flow_cumulants.finalize(numerator,denominator,four=method=="cumulant_4")
    return ratio_estimate(numerator,denominator)


def diagnostics(obs,numerator,denominator):
    ratios=ratio_estimate(numerator,denominator)
    method=obs.flow["method"]
    if method=="event_plane":
        from .flow_resolution import full_resolution
        return {"subevent_correlation":float(ratios[-1]),"resolution":float(full_resolution(ratios[-1]))}
    if method.startswith("cumulant"):
        result={"c2":float(ratios[0])}
        if method=="cumulant_4":
            result["c4"]=float(ratios[1]-2*ratios[0]**2)
        return result
    return {}
