"""Translate the short user configuration into explicit analysis settings."""
import copy
import re
import numpy as np
from pathlib import Path
from .io.config import load_config
from .core.selection import all_charged_alias, normalize_selection

# The left side is what users write in the short configuration.
KINDS = {"pt_spectra": "pt", "mt_spectra": "mt", "yields": "yield"}
CUTS = {"pdg": "pdg", "charged": "charged", "y_range": "y", "eta_range": "eta", "pt_range": "pt",
        "frame": "frame", "input_frame": "input_frame", "cm_rapidity": "cm_rapidity",
        "particle_definition": "particle_definition", "boundary": "boundary"}
BIN_KEYS = {"y_bins": "y", "eta_bins": "eta", "pt_bins": "pt", "mt_bins": "mt"}


def expand_observables(definitions):
    """One short block becomes one calculation per requested particle species.

    The explicit list format from 0.1 remains supported with its original meaning.
    """
    if isinstance(definitions, list):
        return copy.deepcopy(definitions)
    if not isinstance(definitions, dict):
        raise ValueError("observables must be a mapping of modules or an explicit list")
    result = []
    for module_name, settings in definitions.items():
        if not isinstance(settings, dict):
            raise ValueError(f"Settings for {module_name} must be a mapping")
        options = copy.deepcopy(settings)
        enabled = options.pop("enabled", True)
        if not isinstance(enabled, bool):
            raise ValueError(f"{module_name}.enabled must be true or false")
        if not enabled:
            continue
        kind = options.pop("kind", KINDS.get(module_name, module_name))
        selection = options.pop("selection", {})
        if not isinstance(selection, dict):
            raise ValueError("selection must be a mapping")
        for short, full in CUTS.items():
            if short in options:
                if full in selection:
                    raise ValueError(f"Duplicate {full} selection in {module_name}")
                selection[full] = options.pop(short)
        axis = options.get("axis", "y")
        expected_axis = axis if kind in ("v1", "v2", "vn") else {
            "dndy": "y", "net_proton": "y", "dndeta": "eta", "pt": "pt", "invariant_pt": "pt",
            "mt": "mt", "mt_minus_m0": "mt", "invariant_mt": "mt"}.get(kind)
        for bin_key, bin_axis in BIN_KEYS.items():
            if bin_key in options:
                if bin_axis != expected_axis or "bins" in options:
                    raise ValueError(f"Wrong or duplicate bin axis: {module_name}.{bin_key}")
                options["bins"] = options.pop(bin_key)
        if "invariant" in options:
            invariant = options.pop("invariant")
            if not isinstance(invariant, bool) or kind not in ("pt", "mt"):
                raise ValueError("invariant is a boolean option for pt_spectra or mt_spectra")
            if invariant:
                kind = "invariant_" + kind
        if isinstance(options.get("bins"), dict):
            binning = options["bins"]
            if set(binning) != {"min", "max", "count"} or type(binning["count"]) is not int or binning["count"] < 1:
                raise ValueError("Uniform bins require min, max and a positive integer count")
            options["bins"] = np.linspace(binning["min"],binning["max"],binning["count"]+1).tolist()
        name = options.pop("name", module_name)
        definition = dict(options, name=name, kind=kind, selection=selection)
        if all_charged_alias(selection.get('pdg')) and kind != 'ratio':
            definition['selection'] = normalize_selection(selection)
            definition['name'] = f'{name}_all_charged'
            result.append(definition)
            continue
        # Ratios retain their explicit combined numerator/denominator selections.
        if "pdg" in selection and kind != "ratio":
            species = selection["pdg"]
            if not isinstance(species, list) or not species or any(type(pid) is not int for pid in species):
                raise ValueError(f"{module_name}.pdg must be a nonempty list of integer PDG IDs")
            for pdg in species:
                item = copy.deepcopy(definition)
                item["name"] = f"{name}_{pdg}"
                item["selection"]["pdg"] = [pdg]
                result.append(item)
        else:
            result.append(definition)
    return result


def plot_rc_params(settings):
    """Keep figure styling separate from numerical settings and cache identity."""
    if not isinstance(settings, dict) or set(settings)-{"rc_params"}:
        raise ValueError("plot must be a mapping with an optional rc_params mapping")
    params=settings.get("rc_params",{})
    if not isinstance(params,dict) or any(not isinstance(key,str) for key in params):
        raise ValueError("plot.rc_params must map Matplotlib parameter names to values")
    return copy.deepcopy(params)


def load_plot_style(config_path):
    document=load_config(config_path)
    return plot_rc_params(document.get("plot",{}))


def load_job(config_path):
    """Resolve paths relative to the config file; do not change working directory."""
    path = Path(config_path).resolve()
    document = load_config(path)
    if not isinstance(document, dict):
        raise ValueError("Configuration must be a mapping")
    options = copy.deepcopy(document)
    plot_rc_params(options.pop("plot", {}))
    if "input" not in options:
        raise ValueError("Configuration requires input: /path/to/simulation/runs")

    def local_path(value):
        if not isinstance(value, str) or not value:
            raise ValueError("Paths must be nonempty strings")
        return (path.parent / Path(value).expanduser()).resolve()

    source = local_path(options.pop("input"))
    output_name = options.pop("output", "results")
    workers = options.pop("workers", 1)
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    run_pattern = options.pop("run_pattern", None)
    if run_pattern is not None and (not isinstance(run_pattern, str) or not run_pattern):
        raise ValueError("run_pattern must be a nonempty glob pattern")
    if isinstance(output_name, str) and "{run_pattern}" in output_name:
        if run_pattern is None:
            raise ValueError("output with {run_pattern} requires run_pattern")
        # Keep the directory label readable; glob syntax must not become a path.
        label = re.sub(r"[^\w.-]+", "_", run_pattern).strip("_.")
        if not label:
            raise ValueError("run_pattern must contain a name for output with {run_pattern}")
        output_name = output_name.replace("{run_pattern}", label)
    output = local_path(output_name)
    manifest = options.pop("manifest", None)
    if isinstance(manifest, str):
        manifest = load_config(local_path(manifest))
    if manifest is not None and not isinstance(manifest, dict):
        raise ValueError("manifest must be a path or a mapping")
    sampling = options.pop("sampling", None)
    if sampling is not None:
        if not isinstance(sampling, dict):
            raise ValueError("sampling must be a mapping")
        manifest = copy.deepcopy(manifest) if manifest else {}
        manifest.setdefault("defaults", {})["sampling"] = sampling
    options["observables"] = expand_observables(options.get("observables"))
    from .analysis import validate_config
    try:
        validate_config(options)
    except TypeError as exc:
        raise ValueError(f"Invalid observable setting: {exc}") from exc
    return {"root": source, "output": output, "workers": workers, "manifest": manifest, "config": options, "run_pattern": run_pattern}
