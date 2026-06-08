"""Standalone launcher for the chromatoPy desktop wrapper."""

import os
import sys


def _ensure_standard_streams() -> None:
    """PyInstaller windowed Windows apps can start with stdout/stderr set to None."""
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w", encoding="utf-8")


_ensure_standard_streams()

from chromatopy.gui.app import main


if __name__ == "__main__":
    raise SystemExit(main())
