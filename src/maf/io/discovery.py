"""Deterministic recursive discovery without following directory symlinks."""
from dataclasses import dataclass
from pathlib import Path
import logging
import fnmatch
import os
from .metadata import metadata_for


@dataclass
class Run:
    path: Path
    oscar: Path
    metadata: dict


def discover(root, *, rules=(), manifest=None, config_fields=(), recursive=True, pattern=None):
    if pattern is not None and (not isinstance(pattern, str) or not pattern):
        raise ValueError("run_pattern must be a nonempty glob pattern")
    root = Path(root).resolve()
    if not root.is_dir():
        raise ValueError(f"Simulation directory does not exist: {root}")
    found = []
    def onerror(error):
        logging.warning("Discovery could not inspect %s", error)
    for directory, dirs, files in os.walk(root, onerror=onerror, followlinks=False):
        dirs[:] = sorted(d for d in dirs if not d.startswith(".") and d not in ("venv", "__pycache__"))
        name = next((n for n in ("particle_lists.oscar", "particle_lists.oscar.unfinished") if n in files), None)
        if name and (pattern is None or fnmatch.fnmatchcase(Path(directory).name, pattern)):
            path = Path(directory)
            relative = path.relative_to(root).as_posix()
            found.append(Run(path, path/name, metadata_for(path, relative, rules=rules, manifest=manifest, config_fields=config_fields)))
        if not recursive:
            # Include root and immediate children.
            if Path(directory) != root:
                dirs.clear()
    return sorted(found, key=lambda run: str(run.path))
