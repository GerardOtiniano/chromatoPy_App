"""Regression tests for keeping peak selections and their plot artists in sync."""
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from chromatopy.chromatoPy_base import GDGTAnalyzer


class GdgtPeakSelectionTests(unittest.TestCase):
    def setUp(self):
        x = np.linspace(28, 32, 801)
        # Two detected peaks avoid the separate single-peak valley lookup bug.
        y = (1000 * np.exp(-0.5 * ((x - 30) / 0.06) ** 2)
             + 500 * np.exp(-0.5 * ((x - 31) / 0.06) ** 2))
        self.analyzer = GDGTAnalyzer(
            pd.DataFrame({"RT (min)": x, "1022": y}),
            ["1022"], [28, 32], {"1022": "Ia"},
            4000, "synthetic", True, 10, 9, 3, 0.01, 2, None,
        )
        self.analyzer.fig, self.analyzer.axs = self.analyzer.plot_data()
        self.analyzer.current_ax_idx = 0
        self.ax = self.analyzer.axs[0]
        self.messages = []
        self.analyzer.message_callback = self.messages.append

    def tearDown(self):
        plt.close(self.analyzer.fig)

    def click(self, rt):
        self.analyzer.on_click(SimpleNamespace(inaxes=self.ax, xdata=rt))

    def test_repeated_manual_selection_preserves_one_fit_and_area(self):
        self.click(30)
        original = next(iter(self.analyzer.integrated_peaks.values()))
        self.assertGreater(original["area"], 0)
        with patch.object(self.analyzer, "fit_gaussians") as fit:
            self.click(30)
        fit.assert_not_called()
        self.assertEqual(len(self.ax.collections), 1)
        self.assertEqual(len([text for text in self.ax.texts if text.get_gid() != "peak-selection-count"]), 1)
        self.assertEqual(len(self.analyzer.action_stack), 1)
        self.assertIs(next(iter(self.analyzer.integrated_peaks.values())), original)
        self.assertEqual(len(self.analyzer.peak_results["1022"]["Area"]), 1)
        self.analyzer.collect_peak_data()
        self.assertEqual(self.analyzer.peak_results["Ia"]["Area"], [original["area"]])

    def test_undo_after_reselection_removes_both_data_and_shading(self):
        self.click(30)
        self.click(30)
        self.analyzer.undo_last_action()
        self.assertFalse(self.analyzer.integrated_peaks)
        self.assertEqual(len(self.ax.collections), 0)
        self.assertEqual(len([text for text in self.ax.texts if text.get_gid() != "peak-selection-count"]), 0)
        self.assertFalse(self.analyzer.action_stack)
        self.analyzer.collect_peak_data()
        self.assertNotIn("Ia", self.analyzer.peak_results)

    def test_clicking_automatic_selection_does_not_create_an_undo_action(self):
        self.analyzer.reference_peaks = {"Ia": {"Retention Time": [30.0]}}
        self.analyzer.auto_select_peaks()
        original = next(iter(self.analyzer.integrated_peaks.values()))
        self.click(30)
        self.assertFalse(self.analyzer.action_stack)
        self.analyzer.undo_last_action()
        self.assertIs(next(iter(self.analyzer.integrated_peaks.values())), original)
        self.assertEqual(len(self.ax.collections), 1)
        self.assertEqual(len([text for text in self.ax.texts if text.get_gid() != "peak-selection-count"]), 1)
        self.analyzer.collect_peak_data()
        self.assertEqual(self.analyzer.peak_results["Ia"]["Area"], [original["area"]])

    def test_repeated_automatic_selection_does_not_duplicate_shading(self):
        self.analyzer.reference_peaks = {"Ia": {"Retention Time": [30.0]}}
        self.analyzer.auto_select_peaks()
        self.analyzer.auto_select_peaks()
        self.assertEqual(len(self.analyzer.integrated_peaks), 1)
        self.assertEqual(len(self.ax.collections), 1)
        self.assertEqual(len([text for text in self.ax.texts if text.get_gid() != "peak-selection-count"]), 1)
        self.assertEqual(len(self.analyzer.peak_results["1022"]["Area"]), 1)

    def test_undo_after_duplicate_click_still_undoes_latest_real_selection(self):
        self.click(30)
        original = next(iter(self.analyzer.integrated_peaks.values()))
        self.click(31)
        self.click(30)
        self.analyzer.undo_last_action()
        self.assertEqual(list(self.analyzer.integrated_peaks.values()), [original])
        self.assertEqual(len(self.ax.collections), 1)
        self.assertEqual(len([text for text in self.ax.texts if text.get_gid() != "peak-selection-count"]), 1)
        self.assertEqual(len(self.analyzer.action_stack), 1)

    def test_failed_fit_does_not_create_an_undo_action(self):
        with patch.object(self.analyzer, "fit_gaussians", side_effect=RuntimeError("test failure")):
            self.click(30)
        self.assertFalse(self.analyzer.integrated_peaks)
        self.assertFalse(self.analyzer.action_stack)
        self.assertEqual(len(self.ax.collections), 0)
        self.assertTrue(any("test failure" in message for message in self.messages))

    def test_rejected_peak_has_only_its_placeholder_undo_action(self):
        self.analyzer.baseline_threshold["1022"] = 1e9
        self.click(30)
        self.assertEqual(len(self.analyzer.action_stack), 1)
        self.assertEqual(self.analyzer.action_stack[0][0], "add_nopeak")
        self.analyzer.undo_last_action()
        self.assertFalse(self.analyzer.integrated_peaks)
        self.assertFalse(self.analyzer.no_peak_lines)
        self.assertEqual(len([text for text in self.ax.texts if text.get_gid() != "peak-selection-count"]), 0)


if __name__ == "__main__":
    unittest.main()
