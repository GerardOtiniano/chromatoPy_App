"""Validate chromatoPy release version consistency."""

from __future__ import annotations

import os
import re
import sys
import tomllib
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _pyproject_version() -> str:
    with (PROJECT_ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["tool"]["poetry"]["version"]


def _readme_version() -> str:
    readme = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
    match = re.search(r"# chromatoPy App \(v\.?([0-9]+(?:\.[0-9]+)*(?:[A-Za-z0-9.\-+]*)?)\)", readme)
    if not match:
        raise ValueError("README.md must include a heading like '# chromatoPy App (v2.3.3)'.")
    return match.group(1)


def _tag_version() -> str | None:
    tag = os.environ.get("GITHUB_REF_NAME") or os.environ.get("VERSION_TAG")
    if not tag:
        return None
    return tag[1:] if tag.startswith("v") else tag


def main() -> int:
    expected = _pyproject_version()
    observed = {
        "README.md": _readme_version(),
    }
    tag = _tag_version()
    if tag:
        observed["release tag"] = tag

    mismatches = [
        f"{source} has {version}, expected {expected}"
        for source, version in observed.items()
        if version != expected
    ]
    if mismatches:
        print("Version mismatch:")
        for mismatch in mismatches:
            print(f"  - {mismatch}")
        return 1

    print(f"Version check passed: {expected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
