"""Measured caption placement and leader lines follow plot changes."""
from types import SimpleNamespace
import unittest
import numpy as np
from matplotlib.text import Text
import test_gdgt_peak_selection as selection_tests


class GdgtCaptionLayoutTests(unittest.TestCase):
    setUp = selection_tests.GdgtPeakSelectionTests.setUp
    tearDown = selection_tests.GdgtPeakSelectionTests.tearDown
    click = selection_tests.GdgtPeakSelectionTests.click

    def crowded_selection(self):
        self.click(30)
        for rt in (30.01, 30.02, 30.03, 30.04):
            self.analyzer._register_no_peak(self.ax, 0, rt, "1022", "grey")

    def assert_separated(self):
        self.analyzer.fig.canvas.draw()
        renderer = self.analyzer.fig.canvas.get_renderer()
        boxes = [Text.get_window_extent(p["text"], renderer)
                 for p in self.analyzer.integrated_peaks.values() if p["text"].get_visible()]
        boxes.append(self.analyzer.selection_count_texts["1022"].get_window_extent(renderer))
        for i, box in enumerate(boxes):
            for other in boxes[i + 1:]:
                self.assertFalse(box.overlaps(other))
        bounds = self.ax.get_window_extent(renderer)
        for box in boxes:
            self.assertGreaterEqual(box.x0, bounds.x0)
            self.assertLessEqual(box.x1, bounds.x1)
            self.assertGreaterEqual(box.y0, bounds.y0)
            self.assertLessEqual(box.y1, bounds.y1)

    def test_close_captions_separate_after_zoom_and_resize(self):
        self.crowded_selection()
        self.assert_separated()
        self.analyzer.apply_window_bounds(29.8, 30.3)
        self.assert_separated()
        self.analyzer.fig.set_size_inches(6, 4)
        self.assert_separated()
        # Toolbar-style changes also trigger layout, not just window controls.
        self.ax.set_xlim(29.9, 30.2)
        self.ax.set_ylim(0, 1600)
        self.assert_separated()

    def test_leader_ends_five_percent_above_fitted_apex(self):
        self.click(30)
        peak = next(iter(self.analyzer.integrated_peaks.values()))
        for limits in ((0, 1200), (0, 2000)):
            self.ax.set_ylim(*limits)
            self.analyzer.fig.canvas.draw()
            caption = peak["text"]
            fit = peak["fit"]
            apex = int(np.argmax(fit["y"]))
            self.assertAlmostEqual(caption.xy[0], fit["x"][apex])
            self.assertAlmostEqual(caption.xy[1], fit["y"][apex] + 0.05 * (limits[1] - limits[0]))
            self.assertEqual(caption.arrowprops["color"], "grey")
            self.assertIsNotNone(caption.arrow_patch)

    def test_out_of_window_caption_hides_and_returns(self):
        self.click(30)
        caption = next(iter(self.analyzer.integrated_peaks.values()))["text"]
        self.ax.set_xlim(30.5, 31.5)
        self.analyzer.fig.canvas.draw()
        self.assertFalse(caption.get_visible())
        self.ax.set_xlim(29.5, 31.5)
        self.analyzer.fig.canvas.draw()
        self.assertTrue(caption.get_visible())

    def test_d_r_t_remove_captions_and_their_connectors(self):
        self.crowded_selection()
        self.analyzer.on_key(SimpleNamespace(key="d"))
        self.assert_separated()
        for key in ("r", "t"):
            self.analyzer.on_key(SimpleNamespace(key=key))
            self.analyzer.fig.canvas.draw()
            self.assertFalse(self.analyzer.integrated_peaks)
            self.assertEqual(len(self.ax.texts), 1)
            self.assertEqual(self.ax.texts[0].get_text(), "0:1")
            self.click(30)
            self.assert_separated()


if __name__ == "__main__":
    unittest.main()
