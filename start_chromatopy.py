"""Launch chromatoPy from a source checkout.

Run with:
    python start_chromatopy.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _ensure_src_on_path() -> None:
    project_root = Path(__file__).resolve().parent
    src_dir = project_root / "src"
    if src_dir.exists():
        sys.path.insert(0, str(src_dir))


def _ensure_standard_streams() -> None:
    if sys.stdout is None:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
    if sys.stderr is None:
        sys.stderr = open(os.devnull, "w", encoding="utf-8")


def main() -> int:
    os.environ.setdefault("QT_API", "pyside6")
    _ensure_standard_streams()
    _ensure_src_on_path()

    try:
        from chromatopy.gui.app import main as run_chromatopy
    except ModuleNotFoundError as exc:
        missing = exc.name or "a required package"
        print(
            f"Cannot start chromatoPy because {missing!r} is not installed.\n"
            "Install the project dependencies first, then run this launcher again:\n"
            "    poetry install\n"
            "    poetry run python start_chromatopy.py",
            file=sys.stderr,
        )
        return 1

    return run_chromatopy()


if __name__ == "__main__":
    raise SystemExit(main())
