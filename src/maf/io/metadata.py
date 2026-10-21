"""Typed metadata from explicit manifests, configurations and regex rules."""
import copy
import json
import re
from .config import load_config


def value_at(metadata, key):
    current = metadata
    for part in key.split("."):
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def put(metadata, key, value):
    parts = key.split(".")
    for part in parts[:-1]:
        metadata = metadata.setdefault(part, {})
    metadata[parts[-1]] = value


def update_nested(target, source):
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            update_nested(target[key], value)
        else:
            target[key] = copy.deepcopy(value)
    return target


def parse_name(name, rules=()):
    metadata = {}
    # Only universally visible identity fields: no model parameter defaults.
    match = re.match(r"^(\d+(?:\.\d+)?)(?:_([a-z][a-z0-9_]*))?_?([A-Z][a-z]?[A-Z][a-z]?)", name)
    if match:
        metadata["energy"] = {"value": float(match[1]), "unit": "GeV"}
        # Preserve a literal directory label without inferring an IC model.
        if match[2]:
            metadata["directory_label"] = match[2].rstrip("_")
        system = match[3]
        # Legacy AuAU capitalization occurs in the supplied examples.
        if len(system) == 3 and len(name) > match.end() and name[match.end()].isupper():
            system += name[match.end()]
        if len(system) % 2 == 0:
            mid = len(system) // 2
            if system[:mid].lower() == system[mid:].lower():
                element = system[:mid].capitalize()
                system = element + "+" + element
        if "+" not in system:
            parts = re.findall(r"[A-Z][a-z]?", system)
            system = "+".join(parts)
        metadata["collision_system"] = system
    run = re.search(r"(?:^|_)run_*(\d+)(?:_|$)", name)
    if run:
        metadata["run_id"] = int(run[1])
    casters = {"str": str, "int": int, "float": lambda s: float(s.replace("p", ".")),
               "bool": lambda s: {"true": True, "false": False}[s.lower()]}
    for rule in rules:
        match = re.search(rule["pattern"], name)
        if match:
            for group, spec in rule.get("fields", {}).items():
                if match.group(group) is None:
                    continue
                value = casters[spec.get("type", "str")](match.group(group))
                if "unit" in spec:
                    value = {"value": value, "unit": spec["unit"]}
                put(metadata, spec.get("key", group), value)
    return metadata


def metadata_for(run_dir, relative, *, rules=(), manifest=None, config_fields=()):
    result = parse_name(run_dir.name, rules)
    for entry in config_fields:
        path = run_dir / entry["file"]
        if path.exists():
            config = load_config(path)
            for source, target in entry["fields"].items():
                value = value_at(config, source)
                if value is not None:
                    put(result, target, value)
    if manifest:
        update_nested(result, manifest.get("defaults", {}))
        update_nested(result, manifest.get("runs", {}).get(relative, {}))
    result["path"] = str(run_dir.resolve())
    return result


def grouping_key(metadata, keys, identity):
    # An explicit empty list declares one ensemble across all selected files.
    values = [value_at(metadata, key) for key in keys]
    # Missing configurations remain isolated, even from one another.
    if any(v is None for v in values):
        return json.dumps({"unknown_run": identity, "keys": list(keys)}, sort_keys=True)
    return json.dumps(dict(zip(keys, values)), sort_keys=True, separators=(",", ":"))
