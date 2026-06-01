"""Convert IC MS raw text exports into chromatopy-ready CSV files."""

from __future__ import annotations

from pathlib import Path
import csv
import re


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def split_row(line: str):
    line = line.strip()
    if not line:
        return []
    if "\t" in line:
        return [x.strip() for x in line.split("\t")]
    if "," in line:
        return [x.strip() for x in line.split(",")]
    return re.split(r"\s{2,}|\s+", line)


def extract_chromatogram_table(txt_path: Path):
    with txt_path.open("r", encoding="utf-8", errors="ignore") as handle:
        lines = handle.readlines()

    start_idx = None
    for i, line in enumerate(lines):
        if "Chromatogram Data:" in line:
            start_idx = i + 1
            break

    if start_idx is None:
        return None, None

    while start_idx < len(lines) and not lines[start_idx].strip():
        start_idx += 1

    table_rows = []
    for line in lines[start_idx:]:
        stripped = line.strip()
        if not stripped:
            if table_rows:
                break
            continue
        row = split_row(line)
        if not row:
            if table_rows:
                break
            continue
        if len(row) < 2:
            if table_rows:
                break
            continue
        table_rows.append(row)

    if not table_rows:
        return None, None

    first_row = table_rows[0]
    if not all(is_number(x) for x in first_row):
        return first_row, table_rows[1:]
    return None, table_rows


def classify_folder(folder_name: str) -> str:
    lower_name = folder_name.lower()
    if "cation" in lower_name:
        return "cation"
    if "anion" in lower_name:
        return "anion"
    return "uncategorized"


def ic_ms_to_csv(base_path=None, output_base_path=None):
    if base_path is None:
        raise ValueError("A parent directory containing IC MS subfolders is required.")

    parent_dir = Path(base_path)
    if not parent_dir.exists() or not parent_dir.is_dir():
        raise ValueError(f"Provided path is not a valid directory: {parent_dir}")

    output_root = Path(output_base_path) if output_base_path else parent_dir / "Converted Files"
    found_categories = set()
    files_to_process = []

    for subfolder in parent_dir.iterdir():
        if not subfolder.is_dir():
            continue
        if subfolder.name == output_root.name:
            continue
        category = classify_folder(subfolder.name)
        txt_files = list(subfolder.glob("*.txt"))
        if txt_files:
            found_categories.add(category)
            for txt_file in txt_files:
                files_to_process.append((subfolder, category, txt_file))

    if not files_to_process:
        return {
            "input_path": str(parent_dir),
            "output_path": str(output_root),
            "files_exported": [],
            "files_skipped": [],
        }

    output_root.mkdir(exist_ok=True)
    for category in found_categories:
        (output_root / category).mkdir(exist_ok=True)

    exported_files = []
    skipped_files = []
    for subfolder, category, txt_file in files_to_process:
        header, data_rows = extract_chromatogram_table(txt_file)
        if not data_rows:
            skipped_files.append(str(txt_file))
            continue

        output_name = f"{subfolder.name}_{txt_file.stem}.csv"
        output_path = output_root / category / output_name
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            if header:
                writer.writerow(header)
            writer.writerows(data_rows)
        exported_files.append(str(output_path))

    return {
        "input_path": str(parent_dir),
        "output_path": str(output_root),
        "files_exported": exported_files,
        "files_skipped": skipped_files,
    }
