"""Header-driven streaming OSCAR2013/Extended final-particle reader.

The last out block is committed only on its matching event end marker.
Intermediate snapshots and in blocks do not become extra events.
"""
import re
import warnings
from pathlib import Path
import numpy as np
from ..core.event import Event, Particles


class OscarFormatError(ValueError):
    pass


class OscarWarning(UserWarning):
    pass


class OscarReader:
    def __init__(self, path, *, errors="strict", shell_rtol=5e-4, shell_atol=1e-6):
        if errors not in ("strict", "warn"):
            raise ValueError("errors must be strict or warn")
        if not np.isfinite([shell_rtol, shell_atol]).all() or shell_rtol < 0 or shell_atol < 0:
            raise ValueError("Mass-shell tolerances must be finite and nonnegative")
        self.path = Path(path)
        self.errors = errors
        self.shell_rtol, self.shell_atol = shell_rtol, shell_atol
        self.diagnostics = {}
        self.header = {}

    def _issue(self, message, line):
        message = f"{self.path}:{line}: {message}"
        self.diagnostics["issues"] += 1
        if self.errors == "strict":
            raise OscarFormatError(message)
        warnings.warn(message, OscarWarning, stacklevel=3)

    def _particles(self, rows, columns):
        data = np.asarray(rows, dtype=float).reshape((-1, len(columns)))
        if not np.isfinite(data).all():
            raise ValueError("Nonfinite particle record")
        index = {name: i for i, name in enumerate(columns)}
        def col(name):
            return data[:, index[name]]
        for name in ("pdg", "ID", "charge"):
            if name in index and np.any(col(name) != np.floor(col(name))):
                raise ValueError(f"Noninteger {name}")
        required = {"p0", "px", "py", "pz", "mass", "pdg", "ID", "charge"}
        p = Particles(col("p0"), col("px"), col("py"), col("pz"), col("mass"), col("pdg"), col("charge"),
                      col("ID") if "ID" in index else None,
                      {n: col(n) for n in columns if n not in required})
        p.validate_shell(self.shell_rtol, self.shell_atol)
        return p

    def __iter__(self):
        self.diagnostics = {"events": 0, "discarded": 0, "issues": 0, "out_blocks": 0, "unknown_impact": 0}
        columns = None
        key = None
        pending = None
        invalid = False
        seen = set()
        line_number = 0
        # A file, including .unfinished, can contain fully completed events.
        with self.path.open(encoding="utf-8", errors="strict") as stream:
            for line in stream:
                line_number += 1
                tokens = line.split()
                if not tokens:
                    continue
                if tokens[0].startswith("#!"):
                    if columns is not None or tokens[0] not in ("#!OSCAR2013", "#!OSCAR2013Extended") or len(tokens) < 3 or tokens[1] != "particle_lists":
                        raise OscarFormatError(f"{self.path}:{line_number}: unsupported or repeated OSCAR header")
                    columns = tokens[2:]
                    required = {"p0", "px", "py", "pz", "mass", "pdg", "charge"}
                    if len(set(columns)) != len(columns) or not required.issubset(columns):
                        raise OscarFormatError("Missing or duplicated OSCAR column names")
                    self.header = {"format": tokens[0][2:], "columns": columns}
                    continue
                if tokens[:2] != ["#", "event"]:
                    if tokens[0] != "#":
                        self._issue("Particle record outside declared block", line_number)
                        if key is not None:
                            invalid = True
                    elif len(tokens) > 1 and tokens[1] == "Units:":
                        self.header["units"] = tokens[2:]
                        if columns is not None:
                            if len(tokens[2:]) != len(columns):
                                raise OscarFormatError("Units header length does not match columns")
                            for name in ("p0", "px", "py", "pz", "mass"):
                                if tokens[2:][columns.index(name)] != "GeV":
                                    raise OscarFormatError("Only GeV momentum, energy and mass units are supported")
                    elif len(tokens) > 1 and tokens[1].startswith("SMASH"):
                        self.header["generator"] = " ".join(tokens[1:])
                    continue
                if columns is None:
                    raise OscarFormatError("Missing OSCAR2013 column header")
                try:
                    event_id = int(tokens[2])
                    ensemble = int(tokens[tokens.index("ensemble") + 1]) if "ensemble" in tokens else 0
                    kind = next(t for t in tokens[3:] if t in ("in", "out", "end"))
                    current = (event_id, ensemble)
                except (ValueError, IndexError, StopIteration):
                    self._issue("Malformed event marker", line_number)
                    invalid = True
                    continue
                if current != key:
                    if key is not None:
                        self.diagnostics["discarded"] += 1
                        self._issue("Event missing end marker; discarding", line_number)
                    key, pending, invalid = current, None, False
                if kind == "end":
                    try:
                        weight = float(tokens[tokens.index("weight") + 1]) if "weight" in tokens else 1.0
                        impact_raw = float(tokens[tokens.index("impact")+1]) if "impact" in tokens else None
                        # SMASH ModusDefault returns -1 outside collider mode, including List afterburners.
                        # It denotes unavailable collision geometry, not a malformed event or b = 0.
                        impact = None if impact_raw == -1.0 else impact_raw
                        flag = tokens[tokens.index("scattering_projectile_target")+1] if "scattering_projectile_target" in tokens else None
                        if impact is not None and (not np.isfinite(impact) or impact < 0):
                            raise ValueError("Invalid event impact parameter")
                        if flag is not None and flag not in ("yes", "no"):
                            raise ValueError("Invalid scattering_projectile_target flag")
                        if current in seen:
                            raise ValueError("Duplicate event/ensemble identifier")
                        if pending is None:
                            raise ValueError("Event end without an out block")
                        if not invalid:
                            event = Event(pending, event_id, ensemble, weight=weight,
                                          metadata={"end_fields": tokens[3:],
                                              "impact": impact, "impact_raw": impact_raw,
                                              "scattering_projectile_target": flag == "yes" if flag is not None else None})
                            seen.add(current)
                            self.diagnostics["events"] += 1
                            self.diagnostics["unknown_impact"] += int(impact is None)
                            yield event
                        else:
                            self.diagnostics["discarded"] += 1
                    except (ValueError, IndexError) as exc:
                        self.diagnostics["discarded"] += 1
                        self._issue(str(exc), line_number)
                    key, pending, invalid = None, None, False
                    continue
                try:
                    count = int(tokens[tokens.index(kind) + 1])
                    if count < 0:
                        raise ValueError("Negative particle count")
                except (ValueError, IndexError) as exc:
                    self._issue(f"Invalid block count: {exc}", line_number)
                    invalid = True
                    continue
                rows = []
                for _ in range(count):
                    record = stream.readline()
                    line_number += 1
                    if not record or record.lstrip().startswith("#"):
                        # A truncated declared block cannot safely be resynchronized.
                        self.diagnostics["discarded"] += 1
                        self._issue("Truncated particle block; remaining file discarded", line_number)
                        return
                    fields = record.split()
                    try:
                        if len(fields) != len(columns):
                            raise ValueError("Particle column count does not match header")
                        rows.append([float(value) for value in fields])
                    except ValueError as exc:
                        invalid = True
                        self._issue(str(exc), line_number)
                if kind == "out":
                    self.diagnostics["out_blocks"] += 1
                    if not invalid:
                        try:
                            pending = self._particles(rows, columns)
                        except ValueError as exc:
                            invalid = True
                            self._issue(str(exc), line_number)
            if columns is None:
                raise OscarFormatError("Missing OSCAR header")
            if key is not None:
                self.diagnostics["discarded"] += 1
                self._issue("Incomplete final event; discarding", line_number)
