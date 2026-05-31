"""/src/chromatopy/FID/__init__.py"""

__all__ = [
    "integration",
    "peak_label_editor",
    "plot_chromatogram",
    "plot_chromatogram_cluster",
    "cluster",
    "save_results",
    "load_results",
    "delete_samples",
]

name = "FID"


def __getattr__(attribute):
    if attribute == "integration":
        from .FID_integration import integration

        return integration
    if attribute == "peak_label_editor":
        from .peak_labels_editor import peak_label_editor

        return peak_label_editor
    if attribute in {"plot_chromatogram", "plot_chromatogram_cluster"}:
        from .FID_General import plot_chromatogram, plot_chromatogram_cluster

        return {
            "plot_chromatogram": plot_chromatogram,
            "plot_chromatogram_cluster": plot_chromatogram_cluster,
        }[attribute]
    if attribute == "cluster":
        from .bouqueter import cluster

        return cluster
    if attribute in {"save_results", "load_results", "delete_samples"}:
        from .Tools import delete_samples, load_results, save_results

        return {
            "save_results": save_results,
            "load_results": load_results,
            "delete_samples": delete_samples,
        }[attribute]
    raise AttributeError(f"module {__name__!r} has no attribute {attribute!r}")
