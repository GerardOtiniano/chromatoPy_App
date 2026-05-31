"""Public package exports for chromatoPy."""

from __future__ import annotations

from importlib import import_module

__all__ = ["hplc_integration", "hplc_to_csv", "assign_indices", "FID", "IRMS"]

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
