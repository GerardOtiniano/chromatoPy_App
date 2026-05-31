"""Local JSON-backed settings memory for the desktop app."""

from __future__ import annotations

import json
from pathlib import Path


MEMORY_PATH = Path.home() / ".chromatopy_gui_memory.json"


def load_settings_memory() -> dict:
    if not MEMORY_PATH.exists():
        return {"compound_histories": []}
    try:
        return json.loads(MEMORY_PATH.read_text(encoding="utf-8"))
    except Exception:
        return {"compound_histories": []}


def save_settings_memory(memory: dict) -> None:
    MEMORY_PATH.write_text(json.dumps(memory, indent=2), encoding="utf-8")


def list_compound_histories() -> list[list[str]]:
    memory = load_settings_memory()
    histories = memory.get("compound_histories", [])
    return [entry for entry in histories if isinstance(entry, list) and entry]


def remember_compound_history(compounds: list[str]) -> None:
    cleaned = [compound.strip() for compound in compounds if compound.strip()]
    if not cleaned:
        return
    histories = list_compound_histories()
    histories = [entry for entry in histories if entry != cleaned]
    histories.insert(0, cleaned)
    save_settings_memory({"compound_histories": histories[:20]})


def delete_compound_history(compounds: list[str]) -> None:
    cleaned = [compound.strip() for compound in compounds if compound.strip()]
    histories = [entry for entry in list_compound_histories() if entry != cleaned]
    save_settings_memory({"compound_histories": histories})
