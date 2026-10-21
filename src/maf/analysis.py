"""Streaming run analysis, deterministic grouping and source-aware caching."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import logging
import platform
import time
import numpy as np
import warnings
import zipfile
from . import __version__
from .io.oscar import OscarReader
from .io.discovery import discover
from .io.metadata import grouping_key, value_at
from .observables import Observable
from .statistics.aggregation import Result


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def result_directory_name(result):
    return "group-" + digest({"signature": result.signature, "group": result.group})[:16]


def write_result_groups(groups, output, *, failures=(), signature=None):
    """One ensemble lives in results/; several groups get children of results/."""
    output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    index={"schema_version":1,"signature":signature,"groups":[],"summaries":[],"failures":list(failures)}
    for result in groups:
        if failures:
            result.messages.append(f"Analysis is partial: {len(failures)} runs failed; inspect failures.json.")
        directory=output/"results"
        if len(groups)>1:
            directory/=result_directory_name(result)
        result.save(directory/"accumulators.npz")
        result.export(directory)
        index["groups"].append(str(directory/"accumulators.npz"))
        index["summaries"].append((directory/"summary.json").relative_to(output).as_posix())
    patterns = {result.run_pattern for result in groups}
    if len(patterns) == 1:
        index["run_pattern"] = patterns.pop()
    (output/"index.json").write_text(json.dumps(index,indent=2)+"\n")
    return index


def source_fingerprint(path, content_hash=False):
    stat = path.stat()
    info = {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if content_hash:
        hasher = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024*1024), b""):
                hasher.update(block)
        info["sha256"] = hasher.hexdigest()
    return info


def validate_config(config):
    allowed = {"observables", "group_by", "reader", "parsing_rules", "config_fields", "content_hash", "max_events", "event_selection"}
    unknown = set(config) - allowed
    if unknown:
        raise ValueError(f"Unknown configuration keys: {unknown}")
    definitions = config.get("observables", [])
    observables = [Observable(**o) for o in definitions]
    if not observables or len({o.name for o in observables}) != len(observables):
        raise ValueError("Configure at least one observable with unique names")
    keys = config.get("group_by", ["path"])
    if not isinstance(keys, list) or any(not isinstance(k, str) or not k for k in keys):
        raise ValueError("group_by must be a list of metadata keys; [] combines all selected files")
    maximum = config.get("max_events")
    if maximum is not None and (not isinstance(maximum, int) or isinstance(maximum, bool) or maximum < 1):
        raise ValueError("max_events must be a positive integer")
    event_selection = config.get("event_selection", {})
    if "scattering_projectile_target" in event_selection and not isinstance(event_selection["scattering_projectile_target"], bool):
        raise ValueError("scattering_projectile_target must be a boolean")
    # Validate reader options before any worker starts.
    OscarReader("unused", **config.get("reader", {}))
    return observables


def analyze_run(run, config, *, cache=None, resume=True):
    """Analyze one input, returning the same result for cached and fresh work."""
    result, _ = _analyze_run_with_status(run, config, cache=cache, resume=resume)
    return result


def _analyze_run_with_status(run, config, *, cache=None, resume=True):
    observables = validate_config(config)
    signature = digest({"config": config, "version": __version__})
    source_id = digest(str(run.oscar.resolve()))
    before = source_fingerprint(run.oscar, config.get("content_hash", False))
    group = grouping_key(run.metadata, config.get("group_by", ["path"]), source_id)
    cache_key = digest({"source": before, "metadata": run.metadata, "signature": signature})
    cache_path = Path(cache)/f"{cache_key}.npz" if cache else None
    if resume and cache_path and cache_path.exists():
        try:
            result = Result.load(cache_path)
            if result.signature == signature and result.sources[source_id]["fingerprint"] == before:
                logging.debug("Using cache for %s", run.path)
                return result, True
        except (ValueError, KeyError, OSError, EOFError, zipfile.BadZipFile):
            logging.warning("Ignoring invalid cache %s", cache_path)
    result = Result(observables, signature, group, config=config)
    sampling = run.metadata.get("sampling", {})
    mode = sampling.get("mode", "run")
    if mode not in ("event", "hydro", "run"):
        raise ValueError("sampling.mode must be event, hydro or run")
    if mode == "run":
        message = "Hydro independence is unspecified: using runs as independent clusters; event-weighted central values. Independence between runs still requires verification."
        warnings.warn(message, UserWarning)
        result.messages.append(message)
    namespace = str(sampling.get("namespace", source_id))
    hydro_ids = sampling.get("hydro_ids", {})
    reader = OscarReader(run.oscar, **config.get("reader", {"errors": "warn"}))
    max_events = config.get("max_events")
    if max_events is not None and (not isinstance(max_events, int) or max_events < 1):
        raise ValueError("max_events must be a positive integer")
    event_selection = config.get("event_selection", {})
    unknown = set(event_selection) - {"impact", "scattering_projectile_target"}
    if unknown:
        raise ValueError(f"Unknown event selection keys: {unknown}")
    if "impact" in event_selection:
        bounds = event_selection["impact"]
        if len(bounds) != 2 or not (0 <= bounds[0] < bounds[1]):
            raise ValueError("Impact parameter selection requires 0 <= low < high in fm")
    selected_events = 0
    for index,event in enumerate(reader):
        if event_selection:
            include = True
            for key, bounds in event_selection.items():
                value = event.metadata.get(key)
                if value is None:
                    raise ValueError(f"Missing event field required for selection: {key}")
                include &= (bounds[0] <= value < bounds[1]) if key == "impact" else value == bounds
            if not include:
                if max_events is not None and index+1 >= max_events:
                    break
                continue
        selected_events += 1
        identity = f"{event.event_id}:{event.ensemble_id}"
        event_info = sampling.get("events", {}).get(identity, {})
        extra_metadata = event_info.get("metadata", {})
        if not isinstance(extra_metadata, dict):
            raise ValueError("sampling.events metadata must be a mapping")
        if set(extra_metadata) & set(event.metadata):
            raise ValueError("Manifest event metadata cannot overwrite OSCAR event fields")
        event.metadata.update(extra_metadata)
        if "weight" in event_info:
            event.weight = float(event_info["weight"])
            event.__post_init__()
        if mode == "hydro":
            hydro_id = event_info.get("hydro_id", hydro_ids.get(identity, event.hydro_id))
            if hydro_id is None:
                raise ValueError(f"Missing hydro ID for {identity}; cannot estimate independent hydro statistics")
            key = json.dumps(["hydro", namespace, str(hydro_id)])
        else:
            key = json.dumps([mode, source_id, identity if mode == "event" else "run"])
        result.add(key, mode, event, [obs.evaluate(event) for obs in observables])
        if max_events is not None and index+1 >= max_events:
            break
    after = source_fingerprint(run.oscar, config.get("content_hash", False))
    if before != after:
        raise RuntimeError(f"Input changed during analysis: {run.oscar}; no checkpoint written")
    if not result.units:
        raise ValueError(f"No complete usable events in {run.oscar}")
    if reader.diagnostics["unknown_impact"]:
        result.messages.append("Some input events have no impact parameter (SMASH -1 or missing field); no impact cut was applied. See source diagnostics: unknown_impact.")
    if reader.diagnostics["discarded"] or reader.diagnostics["issues"]:
        result.messages.append("Malformed or incomplete input was excluded; inspect source diagnostics before using results.")
    result.sources[source_id] = {"fingerprint": before, "metadata": run.metadata,
        "reader_header": reader.header, "diagnostics": reader.diagnostics,
        "analysis_software": {"package": __version__, "python": platform.python_version(), "numpy": np.__version__},
        "limited": max_events is not None, "selected_events": selected_events, "sampling_mode": mode}
    if cache_path:
        result.save(cache_path)
    return result, False


def merge_results(results):
    groups = {}
    # Source order, independent of process completion order, controls summation.
    for result in sorted(results, key=lambda r: sorted(r.sources)):
        key = (result.signature, result.group)
        if key in groups:
            groups[key].merge(result)
        else:
            groups[key] = result
    return [groups[key] for key in sorted(groups)]


def analyze(root, config, output, *, workers=1, resume=True, dry_run=False, manifest=None, fail_fast=False, run_pattern=None):
    started = time.monotonic()
    logging.info("MAF %s | discovering inputs in %s | pattern=%s", __version__, root, run_pattern or "(all)")
    # Validate in the parent, so worker errors do not obscure configuration errors.
    validate_config(config)
    if workers < 1:
        raise ValueError("workers must be positive")
    runs = discover(root, rules=config.get("parsing_rules", ()), manifest=manifest,
                    config_fields=config.get("config_fields", ()), pattern=run_pattern)
    if not runs:
        raise ValueError("No OSCAR simulation runs found")
    logging.info("Found %d files | workers=%d | resume=%s", len(runs), workers, "on" if resume else "off")
    if dry_run:
        logging.info("Dry run completed | elapsed=%.1f s", time.monotonic() - started)
        return runs
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    results, failures = [], []
    analyzed = cached = 0

    def report():
        completed = analyzed + cached + len(failures)
        return (completed, len(runs), 100 * completed / len(runs), analyzed, cached,
                len(failures), time.monotonic() - started)

    def accept(run, callback):
        nonlocal analyzed, cached
        try:
            result, used_cache = callback()
            results.append(result)
            cached += int(used_cache)
            analyzed += int(not used_cache)
        except Exception as exc:
            logging.error("Failed %s: %s", run.path, exc)
            failures.append({"path": str(run.path), "error": str(exc)})
            if fail_fast:
                logging.info("Stopped: %d/%d files (%.1f%%) | analyzed=%d cached=%d failed=%d | elapsed=%.1f s", *report())
                raise
        # The parent reports both fresh and cached work, including in parallel mode.
        logging.info("Progress: %d/%d files (%.1f%%) | analyzed=%d cached=%d failed=%d | elapsed=%.1f s | %s",
                     *report(), run.path.name)

    if workers == 1:
        for run in runs:
            accept(run, lambda run=run: _analyze_run_with_status(run, config, cache=output/"cache", resume=resume))
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_analyze_run_with_status, run, config, cache=output/"cache", resume=resume): run for run in runs}
            for future in as_completed(futures):
                accept(futures[future], future.result)
    failures.sort(key=lambda failure: failure["path"])
    (output/"failures.json").write_text(json.dumps(failures, indent=2)+"\n")
    if not results:
        raise ValueError("All runs failed; inspect failures.json")
    logging.info("Aggregating %d successful files", len(results))
    groups = merge_results(results)
    for result in groups:
        result.run_pattern = run_pattern
    logging.info("Writing %d result groups to %s", len(groups), output)
    write_result_groups(groups,output,failures=failures,signature=digest({"config":config,"version":__version__}))
    logging.info("Completed: %d/%d files (%.1f%%) | analyzed=%d cached=%d failed=%d | elapsed=%.1f s", *report())
    return groups
