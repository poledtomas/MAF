"""Mergeable sufficient statistics per independent sampling unit."""
from dataclasses import asdict, dataclass, field
from pathlib import Path
import json
import os
import platform
import tempfile
import warnings
import numpy as np
from .errors import jackknife_derived
from .. import __version__
from ..observables import Observable
from ..io.dat import observable_axes, write_dat


def json_clean(value):
    if isinstance(value, dict):
        return {str(k): json_clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def observable_config(obs):
    return {k:v for k,v in asdict(obs).items() if not k.startswith("_")}


@dataclass
class Result:
    observables: list
    signature: str
    group: str
    units: dict = field(default_factory=dict)
    sources: dict = field(default_factory=dict)
    messages: list = field(default_factory=list)
    config: dict = field(default_factory=dict)
    run_pattern: object = None

    @property
    def size(self):
        return sum(o.accumulator_size for o in self.observables)

    def add(self, key, mode, event, values):
        if key not in self.units:
            self.units[key] = {"mode": mode, "a": np.zeros(self.size), "b": np.zeros(self.size), "w": 0.0, "events": 0}
        unit = self.units[key]
        if unit["mode"] != mode:
            raise ValueError("Inconsistent sampling modes for the same unit")
        a = np.concatenate([value[0] for value in values])
        b = np.concatenate([value[1] for value in values])
        unit["a"] += event.weight * a
        unit["b"] += event.weight * b
        unit["w"] += event.weight
        unit["events"] += 1

    def merge(self, other):
        if self.signature != other.signature or self.group != other.group:
            raise ValueError("Cannot merge different analysis definitions or metadata groups")
        if {u["mode"] for u in self.units.values()} != {u["mode"] for u in other.units.values()}:
            raise ValueError("Cannot merge different sampling modes; provide a consistent sampling hierarchy")
        if self.sources.keys() & other.sources.keys():
            raise ValueError("Duplicate source files: refusing to double count")
        if self.run_pattern != other.run_pattern:
            self.run_pattern = None
        for key in sorted(other.units):
            unit = other.units[key]
            if key not in self.units:
                self.units[key] = {k: v.copy() if isinstance(v, np.ndarray) else v for k,v in unit.items()}
            else:
                target = self.units[key]
                if target["mode"] != unit["mode"]:
                    raise ValueError("Incompatible sampling modes")
                for name in ("a", "b", "w", "events"):
                    target[name] += unit[name]
        self.sources.update(other.sources)
        self.messages = sorted(set(self.messages + other.messages))
        return self

    def arrays(self):
        a, b = [], []
        for key in sorted(self.units):
            unit = self.units[key]
            # Average afterburner samples within each hydro event before the ensemble mean.
            scale = unit["w"] if unit["mode"] == "hydro" else 1.0
            a.append(unit["a"] / scale)
            b.append(unit["b"] / scale)
        return np.asarray(a).reshape(-1, self.size), np.asarray(b).reshape(-1, self.size)

    def summary(self):
        a, b = self.arrays()
        messages = list(self.messages)
        if len(self.units) < 2:
            messages.append("Fewer than two independent units: statistical uncertainties are undefined (null).")
        output, offset = {}, 0
        for obs in self.observables:
            end = offset + obs.accumulator_size
            value, error, covariance = jackknife_derived(a[:,offset:end], b[:,offset:end], obs.estimate)
            scale = obs.normalization()
            output[obs.name] = {"definition": observable_config(obs), "value": value*scale,
                "stat_error": error*np.abs(scale), "covariance": covariance*np.outer(scale, scale),
                "normalization": scale, "bins": obs.bins, "units": obs.units,
                "axes": observable_axes(observable_config(obs),obs.units),
                "auxiliary": obs.diagnostics(a[:,offset:end].sum(axis=0),b[:,offset:end].sum(axis=0))}
            offset = end
        return json_clean({"schema_version": 1, "analysis_signature": self.signature,
            "group": json.loads(self.group), "events": sum(u["events"] for u in self.units.values()),
            "independent_units": len(self.units), "sampling_modes": sorted({u["mode"] for u in self.units.values()}),
            "messages": sorted(set(messages)), "sources": self.sources,
            "software": {"package": __version__, "python": platform.python_version(), "numpy": np.__version__},
            "config": self.config, "run_pattern": self.run_pattern, "observables": output})

    def save(self, path):
        """Atomic, pickle-free NPZ checkpoint containing full metadata and raw units."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        keys = sorted(self.units)
        metadata = {"schema_version": 1, "version": __version__, "signature": self.signature,
            "group": self.group, "observables": [observable_config(o) for o in self.observables],
            "sources": self.sources, "messages": self.messages, "config": self.config, "run_pattern": self.run_pattern,
            "feature_layout": [{"name": o.name, "features": o.accumulator_labels} for o in self.observables],
            "unit_keys": keys, "modes": [self.units[k]["mode"] for k in keys]}
        fd, temporary = tempfile.mkstemp(prefix=".analysis-", suffix=".npz", dir=path.parent)
        try:
            with os.fdopen(fd, "wb") as stream:
                np.savez_compressed(stream, metadata=json.dumps(json_clean(metadata), sort_keys=True, allow_nan=False),
                    a=np.asarray([self.units[k]["a"] for k in keys]).reshape(-1,self.size),
                    b=np.asarray([self.units[k]["b"] for k in keys]).reshape(-1,self.size),
                    w=np.asarray([self.units[k]["w"] for k in keys]),
                    events=np.asarray([self.units[k]["events"] for k in keys]))
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as archive:
            metadata = json.loads(str(archive["metadata"]))
            if metadata.get("schema_version") != 1 or metadata.get("version") not in ("0.1.0", "0.2.0", "0.2.1", "0.3.0", "0.3.1", "0.3.2", "0.4.0", "0.4.1", __version__):
                raise ValueError("Unsupported result schema or software version; recompute source results")
            result = cls([Observable(**o) for o in metadata["observables"]], metadata["signature"], metadata["group"],
                         sources=metadata["sources"], messages=metadata["messages"], config=metadata["config"],
                         run_pattern=metadata.get("run_pattern"))
            for i,key in enumerate(metadata["unit_keys"]):
                result.units[key] = {"mode": metadata["modes"][i], "a": archive["a"][i].copy(),
                    "b": archive["b"][i].copy(), "w": float(archive["w"][i]), "events": int(archive["events"][i])}
        return result

    def export(self, directory):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        summary = self.summary()
        (directory/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
        for name, result in summary["observables"].items():
            write_dat(directory/f"{name}.dat",name,result,signature=summary["analysis_signature"],run_pattern=self.run_pattern)
            # Replace the package's previous CSV export after writing its DAT replacement.
            (directory/f"{name}.csv").unlink(missing_ok=True)
