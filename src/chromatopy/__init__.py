"""Public package exports for chromatoPy."""

from __future__ import annotations

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import tomllib


def _source_tree_version() -> str:
    pyproject_path = Path(__file__).resolve().parents[2] / "pyproject.toml"
    with pyproject_path.open("rb") as handle:
        return tomllib.load(handle)["tool"]["poetry"]["version"]

try:
    __version__ = version("chromatopy")
except PackageNotFoundError:  # pragma: no cover - source tree without installed metadata
    try:
        __version__ = _source_tree_version()
    except Exception:
        __version__ = "unknown"

__all__ = ["hplc_integration", "hplc_to_csv", "assign_indices", "FID", "IRMS", "__version__"]

name = "chromatoPy"


def __getattr__(name: str):
    if name in {"hplc_integration", "hplc_to_csv", "assign_indices"}:
        module = import_module(".chromatoPy", __name__)
        return getattr(module, name)
    if name == "FID":
        return import_module(".FID", __name__)
    if name == "IRMS":
        return import_module(".IRMS", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | {"FID", "IRMS"})
