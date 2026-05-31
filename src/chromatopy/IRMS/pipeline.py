"""Batch-friendly IRMS processing wrappers."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

from .dxf_reader import extract_all_raw_datasets


def _save_preview_plot(all_df, preview_path: Path) -> None:
    value_columns = [col for col in all_df.columns if col.startswith("v")]
    if not value_columns:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    x_values = all_df["time_s"] / 60.0
    for column in value_columns:
        ax.plot(x_values, all_df[column], linewidth=1.0, label=column)

    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Signal")
    ax.set_title("IRMS Raw Voltage Preview")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(preview_path, dpi=200)
    plt.close(fig)


def process_irms_file(
    file_path: str | Path,
    output_dir: str | Path,
    export_individual_datasets: bool = True,
    export_preview_plot: bool = True,
):
    """Decode a Thermo Isodat file and export batch-friendly outputs."""

    file_path = Path(file_path).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    file_output_dir = output_dir / file_path.stem
    file_output_dir.mkdir(parents=True, exist_ok=True)

    datasets, all_df = extract_all_raw_datasets(str(file_path), verbose=False)

    combined_csv = file_output_dir / f"{file_path.stem}_raw_voltage_data.csv"
    all_df.to_csv(combined_csv, index=False)

    dataset_exports = []
    if export_individual_datasets:
        dataset_dir = file_output_dir / "datasets"
        dataset_dir.mkdir(exist_ok=True)
        for dataset_name, dataset_df in datasets.items():
            dataset_path = dataset_dir / f"{dataset_name}.csv"
            dataset_df.to_csv(dataset_path, index=False)
            dataset_exports.append(str(dataset_path))

    preview_plot = None
    if export_preview_plot:
        preview_plot = file_output_dir / f"{file_path.stem}_preview.png"
        _save_preview_plot(all_df, preview_plot)

    manifest = {
        "input_file": str(file_path),
        "output_directory": str(file_output_dir),
        "combined_csv": str(combined_csv),
        "preview_plot": str(preview_plot) if preview_plot else None,
        "dataset_exports": dataset_exports,
        "datasets": list(datasets.keys()),
        "rows_exported": int(len(all_df)),
    }

    manifest_path = file_output_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    manifest["manifest_path"] = str(manifest_path)
    return manifest


def process_irms_files(
    file_paths: list[str] | tuple[str, ...],
    output_dir: str | Path,
    export_individual_datasets: bool = True,
    export_preview_plot: bool = True,
):
    """Process one or more IRMS files into CSV/PNG/JSON outputs."""

    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    processed = []
    failed = []
    for file_path in sorted(str(Path(path).expanduser()) for path in file_paths):
        try:
            processed.append(
                process_irms_file(
                    file_path=file_path,
                    output_dir=output_dir,
                    export_individual_datasets=export_individual_datasets,
                    export_preview_plot=export_preview_plot,
                )
            )
        except Exception as exc:  # pragma: no cover - exercised in the GUI
            failed.append({"input_file": file_path, "error": str(exc)})

    summary = {
        "output_directory": str(output_dir),
        "processed": processed,
        "failed": failed,
        "processed_count": len(processed),
        "failed_count": len(failed),
    }

    summary_path = output_dir / "irms_processing_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    summary["summary_path"] = str(summary_path)
    return summary
