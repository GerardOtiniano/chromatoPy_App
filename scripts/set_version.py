"""Update chromatoPy app version references."""

from __future__ import annotations

import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
VERSION_PATTERN = re.compile(r"^[0-9]+(?:\.[0-9]+)*(?:[A-Za-z0-9.\-+]*)?$")


def _normalize_version(raw_version: str) -> str:
    version = raw_version.strip()
    if version.startswith("v"):
        version = version[1:]
    if not VERSION_PATTERN.fullmatch(version):
        raise ValueError(f"Invalid version: {raw_version!r}")
    return version


def _replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"Could not update {label}.")
    return updated


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 1:
        print("Usage: python scripts/set_version.py 2.3.4")
        return 2

    version = _normalize_version(args[0])

    pyproject_path = PROJECT_ROOT / "pyproject.toml"
    pyproject = pyproject_path.read_text(encoding="utf-8")
    pyproject = _replace_once(
        pyproject,
        r'^version = "[^"]+"$',
        f'version = "{version}"',
        "pyproject.toml",
    )
    pyproject_path.write_text(pyproject, encoding="utf-8")

    readme_path = PROJECT_ROOT / "README.md"
    readme = readme_path.read_text(encoding="utf-8")
    readme = _replace_once(
        readme,
        r"^# chromatoPy App \(v\.?[0-9]+(?:\.[0-9]+)*(?:[A-Za-z0-9.\-+]*)?\)$",
        f"# chromatoPy App (v{version})",
        "README.md",
    )
    readme_path.write_text(readme, encoding="utf-8")

    print(f"Updated version to {version}. Release tag should be v{version}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
