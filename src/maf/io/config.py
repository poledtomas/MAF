"""JSON always available; TOML (Python 3.11+) and YAML optional."""
import json
from pathlib import Path


def load_config(path):
    path = Path(path)
    with path.open("rb") as stream:
        if path.suffix == ".toml":
            try:
                import tomllib
            except ImportError as exc:
                raise ImportError("TOML requires Python 3.11+; use JSON or optional YAML on older Python") from exc
            return tomllib.load(stream)
        if path.suffix in (".yaml", ".yml"):
            try:
                import yaml
            except ImportError as exc:
                raise ImportError("Install maf[yaml] for YAML") from exc
            try:
                return yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                raise ValueError(f"Invalid YAML in {path}: {exc}") from exc
        return json.load(stream)
