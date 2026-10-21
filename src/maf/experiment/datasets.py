"""Local CSV+metadata and HEPData-style YAML without publication defaults."""
from dataclasses import dataclass
from pathlib import Path
import csv
import numpy as np
from ..io.config import load_config
from ..observables import Observable


@dataclass
class Dataset:
    metadata: dict
    low: np.ndarray
    high: np.ndarray
    value: np.ndarray
    stat: object = None
    systematic: object = None

    def __post_init__(self):
        for key in ("experiment", "reference", "collision_system", "energy", "centrality", "definition"):
            if key not in self.metadata or self.metadata[key] is None:
                raise ValueError(f"Experimental metadata requires {key}")
        definition = self.metadata["definition"]
        required = ("name", "kind", "selection", "units")
        if any(key not in definition for key in required):
            raise ValueError("Experimental definition requires name, kind, selection and units")
        for key in ("particle_definition", "frame"):
            if key not in definition["selection"]:
                raise ValueError(f"Experimental definition.selection requires explicit {key}")
        cut = Observable(**definition)._cut
        if cut.frame == "input" and cut.input_frame == "unknown":
            raise ValueError("Experimental acceptance requires a known lab or cm frame")
        for key in ("low", "high", "value"):
            setattr(self, key, np.asarray(getattr(self, key), dtype=float))
        if self.low.ndim != 1 or self.value.shape != self.low.shape or self.high.shape != self.low.shape or not len(self.value):
            raise ValueError("Experimental bin/value shapes do not match")
        if not all(np.isfinite(v).all() for v in (self.low,self.high,self.value)) or np.any(self.high <= self.low):
            raise ValueError("Invalid experimental bin or value")
        if np.any(self.low[1:] < self.high[:-1]):
            raise ValueError("Experimental bins overlap or are unsorted")
        for key in ("stat", "systematic"):
            error = getattr(self, key)
            if error is not None:
                error = np.asarray(error, dtype=float)
                if error.shape == self.value.shape:
                    error = np.column_stack((error,error))
                if error.shape != (len(self.value),2) or not np.isfinite(error).all() or np.any(error < 0):
                    raise ValueError(f"Invalid {key} errors; expected N or Nx2 nonnegative values")
                setattr(self,key,error)


def load_csv(path, metadata):
    meta = load_config(metadata) if isinstance(metadata,(str,Path)) else metadata
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    def column(name):
        return [float(row[name]) for row in rows]
    def error(name):
        if rows and name in rows[0]:
            return column(name)
        if rows and name+"_low" in rows[0] and name+"_high" in rows[0]:
            return np.column_stack((column(name+"_low"), column(name+"_high")))
        return None
    return Dataset(meta, column("low"), column("high"), column("value"), error("stat"), error("systematic"))


def load_hepdata(path, metadata, *, dependent=0):
    """Read one HEPData dependent variable. Supply physics metadata explicitly.

    Only stat/sys labelled errors are accepted; ambiguous labels are rejected.
    Percent uncertainties and asymmetric errors are preserved.
    """
    table = load_config(path)
    meta = load_config(metadata) if isinstance(metadata,(str,Path)) else metadata
    if len(table["independent_variables"]) != 1:
        raise ValueError("Only one-dimensional HEPData tables are supported")
    bins = table["independent_variables"][0]["values"]
    if any("low" not in b or "high" not in b for b in bins):
        raise ValueError("HEPData point-only x values need explicit bin boundaries")
    values = table["dependent_variables"][dependent]["values"]
    parsed = {"stat": [], "systematic": []}
    def magnitude(value, central):
        return abs(float(value[:-1])*central/100) if isinstance(value,str) and value.endswith("%") else abs(float(value))
    for item in values:
        errors = {}
        for entry in item.get("errors", []):
            label = entry.get("label", "").lower()
            key = "stat" if label in ("stat", "statistical") else "systematic" if label in ("sys", "syst", "systematic") else None
            if key is None or key in errors:
                raise ValueError(f"Ambiguous HEPData error label: {label}; prepare an explicit local dataset")
            if "symerror" in entry:
                err = magnitude(entry["symerror"],item["value"])
                errors[key] = [err,err]
            else:
                pair = entry["asymerror"]
                errors[key] = [magnitude(pair["minus"],item["value"]),magnitude(pair["plus"],item["value"])]
        for key in parsed:
            parsed[key].append(errors.get(key))
    for key in parsed:
        if any(v is None for v in parsed[key]) and not all(v is None for v in parsed[key]):
            raise ValueError("Inconsistent HEPData errors across bins")
    return Dataset(meta, [b["low"] for b in bins], [b["high"] for b in bins], [v["value"] for v in values],
                   **{k: None if all(v is None for v in entries) else entries for k,entries in parsed.items()})
