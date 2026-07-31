import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
from scipy.integrate import simpson


MODULE_PATH = (
    Path(__file__).parents[1]
    / "src"
    / "chromatopy"
    / "FID"
    / "FID_Integration_functions.py"
)
SPEC = importlib.util.spec_from_file_location(
    "fid_integration_functions_under_test", MODULE_PATH)
fid = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fid)


class FidAsymmetricBoundaryTests(unittest.TestCase):
    def test_strong_left_skew_reaches_equal_low_intensity_tails(self):
        amplitude = 295.3794429680486
        center = 10.697536128518585
        width = 0.04113523100916537
        alpha = -8.490299175312083
        x_min, x_max = fid.calculate_gaus_extension_limits(
            center, width, factor=2)
        x, y = fid.extrapolate_gaussian(
            np.array([]), amplitude, center, width, alpha,
            x_min, x_max, step=0.001)

        left, right = fid.calculate_asymmetric_boundaries(
            y, relative_height=0.001)
        peak = int(np.argmax(y))
        endpoint_fractions = y[[left, right]] / np.max(y)
        retained_area = simpson(y=y[left:right + 1], x=x[left:right + 1])
        complete_area = simpson(y=y, x=x)

        self.assertLess(left, peak)
        self.assertGreater(right, peak)
        self.assertTrue(np.all(endpoint_fractions <= 0.0011))
        self.assertGreater(retained_area / complete_area, 0.9997)
        self.assertGreater(x[peak] - x[left], x[right] - x[peak])

    def test_invalid_relative_height_is_rejected(self):
        with self.assertRaises(ValueError):
            fid.calculate_asymmetric_boundaries([0.0, 1.0, 0.0], 0.0)

    def test_asymmetric_fit_integrates_complete_curve_before_display_crop(self):
        x = pd.Series(np.linspace(10.3, 11.1, 1601))
        y = pd.Series(fid.skewed_gaussian(
            x.to_numpy(), 295.3794429680486, 10.697536128518585,
            0.04113523100916537, -8.490299175312083))
        peak = int(np.argmax(y.to_numpy()))

        fit_x, fit_y, area, result = fid.fit_gaussians(
            x, y, peak, [peak], [9, 3], 0.01, 4000,
            mode="asymmetric")
        displayed_area = simpson(y=fit_y, x=fit_x)

        self.assertEqual(result["model_type"], "Asymmetric-Gaussian")
        self.assertTrue(np.all(
            np.asarray([fit_y[0], fit_y[-1]]) / np.max(fit_y) <= 0.0011))
        self.assertGreater(area, displayed_area)
        self.assertGreater(displayed_area / area, 0.9997)


if __name__ == "__main__":
    unittest.main()
