"""The peak counter is informational and follows selections and keyboard controls."""
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import matplotlib.pyplot as plt
import test_gdgt_peak_selection as selection_tests


class GdgtPeakCounterTests(unittest.TestCase):
    setUp = selection_tests.GdgtPeakSelectionTests.setUp
    tearDown = selection_tests.GdgtPeakSelectionTests.tearDown
    click = selection_tests.GdgtPeakSelectionTests.click

    def key(self, key):
        self.analyzer.on_key(SimpleNamespace(key=key))

    def assert_count(self, selected, expected=1, trace="1022", ax=None):
        ax = ax if ax is not None else self.ax
        labels = [t for t in ax.texts if t.get_gid() == "peak-selection-count"]
        self.assertEqual(len(labels), 1)
        self.assertIs(labels[0], self.analyzer.selection_count_texts[trace])
        self.assertEqual(labels[0].get_text(), f"{selected}:{expected}")

    def two_traces(self):
        a = self.analyzer
        plt.close(a.fig)
        a.df["1036"] = a.df["1022"]
        a.traces = ["1022", "1036"]
        a.GDGT_dict = {"1022": ["IIIa", "IIIa''", "IIIa'"], "1036": ["IIa", "IIa'"]}
        a.fig, a.axs = a.plot_data()
        self.ax = a.axs[0]

    def test_initial_counter_is_faint_grey_and_anchored_to_axes(self):
        self.assert_count(0)
        text = self.analyzer.selection_count_texts["1022"]
        self.assertEqual(text.get_position(), (0.98, 0.98))
        self.assertIs(text.get_transform(), self.ax.transAxes)
        self.assertEqual(text.get_color(), "grey")
        self.assertLess(text.get_alpha(), 1)
        self.assertEqual(text.get_horizontalalignment(), "right")
        self.analyzer.apply_window_bounds(29, 31)
        self.assertIs(text.get_transform(), self.ax.transAxes)
        self.assert_count(0)

    def test_manual_peaks_and_absences_count_without_replacement(self):
        self.click(29.7)
        self.assert_count(1)
        self.click(30)
        self.assert_count(2)
        self.assertEqual(len(self.analyzer.no_peak_lines), 1)
        self.click(30)
        self.assert_count(2)  # Duplicate fitted selections remain ignored.
        self.key("d")
        self.assert_count(1)
        self.assertEqual(len(self.analyzer.no_peak_lines), 1)
        self.key("d")
        self.assert_count(0)
        self.key("d")
        self.assert_count(0)

    def test_five_entries_show_five_of_three_without_changing_export(self):
        a = self.analyzer
        a.GDGT_dict = {"1022": ["IIIa", "IIIa''", "IIIa'"]}
        for rt in (31.7, 31, 30.5, 30, 29.7):
            self.click(rt)
        self.assert_count(5, 3)
        self.key("enter")
        self.assertTrue(a.finished)
        self.assertEqual([a.peak_results[name]["Retention Time"][0]
                          for name in a.GDGT_dict["1022"]], [29.7, 30, 30.5])
        self.assertEqual(a.peak_results["IIIa"]["Area"], [0.0])
        self.assertGreater(a.peak_results["IIIa''"]["Area"][0], 0)

    def test_automatic_selection_counts_real_and_absent_peaks(self):
        a = self.analyzer
        a.GDGT_dict = {"1022": ["IIIa", "IIIa''", "IIIa'"]}
        a.reference_peaks = {
            "IIIa'": {"Retention Time": [31.]},
            "IIIa''": {"Retention Time": [30.5]},
            "IIIa": {"Retention Time": [30.]},
        }
        a.auto_select_peaks()
        self.assert_count(3, 3)
        self.key("d")  # The automatic absence marker has a deletion-history entry.
        self.assert_count(2, 3)
        self.key("t")
        self.assert_count(0, 3)
        self.assertTrue(a.t_pressed)

    def test_r_clears_only_active_trace_and_preserves_other_history(self):
        self.two_traces()
        self.click(30)
        other = self.analyzer.axs[1]
        self.analyzer.on_click(SimpleNamespace(inaxes=other, xdata=29.7))
        self.assert_count(1, 3)
        self.assert_count(1, 2, "1036", other)
        self.key("r")
        self.assert_count(0, 3)
        self.assert_count(1, 2, "1036", other)
        self.key("d")
        self.assert_count(0, 3)
        self.assert_count(0, 2, "1036", other)
        self.click(31)
        self.assert_count(1, 3)
        self.key("d")
        self.assert_count(0, 3)

    def test_t_resets_every_trace_and_repeated_resets_do_not_duplicate_labels(self):
        self.two_traces()
        other = self.analyzer.axs[1]
        self.click(30)
        self.analyzer.on_click(SimpleNamespace(inaxes=other, xdata=29.7))
        for key in ("t", "d", "r", "t", "r"):
            self.key(key)
            self.assert_count(0, 3)
            self.assert_count(0, 2, "1036", other)
        self.assertFalse(self.analyzer.integrated_peaks)
        self.assertFalse(self.analyzer.action_stack)
        self.assertFalse(self.analyzer.no_peak_lines)
        self.click(30)
        self.assert_count(1, 3)

    def test_failed_integration_does_not_increment_counter(self):
        with patch.object(self.analyzer, "fit_gaussians", side_effect=RuntimeError("test failure")):
            self.click(30)
        self.assert_count(0)

    def test_below_threshold_placeholder_counts_once(self):
        self.analyzer.baseline_threshold["1022"] = 1e9
        self.click(30)
        self.assert_count(1)
        self.assertEqual(len(self.analyzer.no_peak_lines), 1)
        self.key("d")
        self.assert_count(0)

    def test_existing_single_channel_count_validation_is_preserved(self):
        self.analyzer.schema_type = "single_channel"
        self.key("enter")
        self.assertFalse(self.analyzer.finished)
        self.assert_count(0)


if __name__ == "__main__":
    unittest.main()
