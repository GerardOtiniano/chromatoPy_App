"""Schema detection and configuration helpers for chromatographic inputs."""

from __future__ import annotations

from dataclasses import dataclass
import os
import re

import pandas as pd


TIME_COLUMN_PATTERNS = (
    re.compile(r"^rt\b", re.IGNORECASE),
    re.compile(r"retention", re.IGNORECASE),
    re.compile(r"\btime\b", re.IGNORECASE),
)
EXCLUDED_COLUMNS = {"Sample Name", "rt_corr"}


@dataclass
class DataSchema:
    schema_type: str
    time_column: str
    signal_columns: list[str]
    csv_files: list[str]
    folder_path: str


def find_time_column(columns: list[str]) -> str:
    """Pick the most likely time/retention-time column."""

    for column in columns:
        if column == "RT (min)":
            return column
    for column in columns:
        for pattern in TIME_COLUMN_PATTERNS:
            if pattern.search(column):
                return column
    raise ValueError("Could not identify a time column in the selected data files.")


def list_csv_files(folder_path: str) -> list[str]:
    """Return numerically sorted CSV files from a folder."""

    def sort_key(filename: str):
        numbers = re.findall(r"\d+", filename)
        return int(numbers[0]) if numbers else 0

    return sorted(
        [name for name in os.listdir(folder_path) if name.lower().endswith(".csv")],
        key=sort_key,
    )


def detect_data_schema(folder_path: str) -> DataSchema:
    """Inspect the first CSV and infer whether the data is multi- or single-channel."""

    csv_files = list_csv_files(folder_path)
    if not csv_files:
        raise ValueError("No CSV files were found in the selected folder.")

    preview_path = os.path.join(folder_path, csv_files[0])
    preview_df = pd.read_csv(preview_path, nrows=5)
    cleaned_columns = [col[:-2] if col.endswith(".0") else col for col in preview_df.columns]
    preview_df.columns = cleaned_columns

    time_column = find_time_column(list(preview_df.columns))
    signal_columns = [
        column
        for column in preview_df.columns
        if column not in EXCLUDED_COLUMNS and column != time_column and not str(column).startswith("Unnamed:")
    ]
    if not signal_columns:
        raise ValueError("No signal columns were detected in the selected data files.")

    schema_type = "multi_channel" if len(signal_columns) > 1 else "single_channel"
    return DataSchema(
        schema_type=schema_type,
        time_column=time_column,
        signal_columns=signal_columns,
        csv_files=csv_files,
        folder_path=folder_path,
    )


def build_single_channel_meta(
    signal_column: str,
    compound_names: list[str],
    window_bounds: list[float],
):
    """Build a GDGT-like metadata structure for generic single-channel workflows."""

    cleaned = [name.strip() for name in compound_names if name.strip()]
    if not cleaned:
        cleaned = [signal_column]

    return {
        "names": [["Compounds"]],
        "GDGT_dict": [{signal_column: cleaned if len(cleaned) > 1 else cleaned[0]}],
        "Trace": [[signal_column]],
        "window": [window_bounds],
    }
