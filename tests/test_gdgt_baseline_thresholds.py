"""Detection thresholds belong to a trace within a single sample."""
import contextlib
import importlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from chromatopy.chromatoPy_base import GDGTAnalyzer
from chromatopy.gui.logic import _iter_hplc_peak_entries, calculate_hplc_peak_area_confidence_intervals

workflow = importlib.import_module('chromatopy.hplc_integration')
# Importing the GUI package selects QtAgg; keep these noninteractive checks headless.
plt.switch_backend('Agg')


def data(scale=1):
    x = np.linspace(28, 32, 801)
    y = scale * (1000 * np.exp(-0.5 * ((x - 30) / .06)**2)
                 + 500 * np.exp(-0.5 * ((x - 31) / .06)**2))
    return pd.DataFrame({'RT (min)': x, '1022': y, '1018': y * 1000})


def analyzer(scale=1, traces=None, minimum=None, reference=None):
    traces = traces or ['1022', '1018']
    a = GDGTAnalyzer(data(scale), traces, [28, 32], {'1022': 'Ia', '1018': 'Ic'},
                     4000, 'synthetic', reference is None, 10, 9, 3, .01, 2, None,
                     min_PA=minimum, reference_peaks=reference)
    return a


def select_first_peak(a, index=0):
    a.on_click(SimpleNamespace(inaxes=a.axs[index], xdata=30.0))


class BaselineThresholdTests(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_strong_neighbor_does_not_reject_valid_peak(self):
        a = analyzer()
        self.assertEqual(a.baseline_threshold, {})
        a.fig, a.axs = a.plot_data()
        self.assertLess(a.baseline_threshold['1022'], 1000)
        self.assertGreater(a.baseline_threshold['1018'], 1000)
        select_first_peak(a)
        a.collect_peak_data()
        self.assertGreater(a.peak_results['Ia']['Area'][0], 140)
        self.assertNotEqual(a.peak_results['Ia']['Model Type'], ['Absent'])

    def test_trace_order_does_not_change_threshold_or_area(self):
        areas = []
        thresholds = []
        for traces in (['1022', '1018'], ['1018', '1022']):
            a = analyzer(traces=traces)
            a.fig, a.axs = a.plot_data()
            select_first_peak(a, traces.index('1022'))
            areas.append(next(iter(a.integrated_peaks.values()))['area'])
            thresholds.append(a.baseline_threshold.copy())
        self.assertEqual(thresholds[0], thresholds[1])
        self.assertAlmostEqual(areas[0], areas[1])

    def test_new_sample_computes_own_thresholds_despite_reference_peaks(self):
        first = analyzer(scale=100)
        first.fig, first.axs = first.plot_data()
        previous = first.baseline_threshold.copy()
        second = analyzer(reference={'Ia': {'Retention Time': [30.]}})
        self.assertEqual(second.baseline_threshold, {})
        second.fig, second.axs = second.plot_data()
        self.assertIsNot(first.baseline_threshold, second.baseline_threshold)
        self.assertAlmostEqual(previous['1022'] / second.baseline_threshold['1022'], 100)
        second.auto_select_peaks()
        self.assertGreater(next(iter(second.integrated_peaks.values()))['area'], 140)
        self.assertEqual(first.baseline_threshold, previous)

    def test_r_and_t_recompute_only_the_appropriate_sample_trace(self):
        a = analyzer()
        a.fig, a.axs = a.plot_data()
        original = a.baseline_threshold.copy()
        a.current_ax_idx = 0
        a.baseline_threshold['1022'] = -123
        a.on_key(SimpleNamespace(key='r'))
        self.assertEqual(a.baseline_threshold, original)
        a.baseline_threshold['1018'] = -456
        a.on_key(SimpleNamespace(key='t'))
        self.assertEqual(a.baseline_threshold, original)

    def test_user_amplitude_override_is_recorded_for_each_trace(self):
        a = analyzer(minimum=25)
        a.fig, a.axs = a.plot_data()
        self.assertEqual(a.baseline_threshold, {'1022': 25.0, '1018': 25.0})

    def test_json_metadata_is_not_mistaken_for_a_compound(self):
        peak = {'Area': [148.7], 'Retention Time': [30.]}
        for sample in (
            {'Sample Name': 'old', 'brGDGTs': {'Ia': peak, 'Ic': 0}},
            {'Sample Name': 'new', 'Baseline Thresholds': {'1022': 14., '1018': 0.},
             'brGDGTs': {'Ia': peak, 'Ic': 0}},
        ):
            self.assertEqual(list(_iter_hplc_peak_entries(sample)), [('Ia', peak), ('Ic', None)])

    def test_real_export_saves_separate_thresholds_for_each_sample(self):
        snapshots = []
        references = []

        def noninteractive_run(a):
            # Use actual baseline estimation, fitting, and result collection;
            # replace only the interactive window/wait with deterministic clicks.
            self.assertEqual(a.baseline_threshold, {})
            references.append(a.reference_peaks)
            a.fig, a.axs = a.plot_data()
            for i in range(len(a.traces)):
                select_first_peak(a, i)
            a.collect_peak_data()
            snapshots.append(a.baseline_threshold.copy())
            return a.peak_results, a.fig, a.peak_results, False

        with TemporaryDirectory() as folder:
            root = Path(folder)
            data().to_csv(root/'sample1.csv', index=False)
            data(2).to_csv(root/'sample2.csv', index=False)
            metadata = {'names': [['brGDGTs']], 'Trace': [['1022', '1018']],
                        'window': [[28., 32.]], 'GDGT_dict': [{'1022': 'Ia', '1018': 'Ic'}]}
            with patch.object(GDGTAnalyzer, 'run', noninteractive_run), contextlib.redirect_stdout(io.StringIO()):
                result = workflow.hplc_integration(folder_path=folder, gdgt_meta_set=metadata,
                    edit_metadata=False, schema_type='multi_channel', normalize_by_standard=False)
            samples = [json.loads((Path(result['sample_path'])/f'sample{i}.json').read_text())
                       for i in (1, 2)]
            for sample, thresholds in zip(samples, snapshots):
                self.assertEqual(sample['Baseline Thresholds'], thresholds)
                self.assertGreater(sample['brGDGTs']['Ia']['Area'][0], 140)
            self.assertAlmostEqual(snapshots[1]['1022'] / snapshots[0]['1022'], 2.)
            for reference in references:
                if reference is not None:
                    self.assertNotIn('Baseline Thresholds', reference)
            table = pd.read_csv(result['results_file_path'])
            self.assertEqual(set(table.columns), {'Sample Name', 'Ia', 'Ic'})
            ci = calculate_hplc_peak_area_confidence_intervals(result['output_folder'])
            ci_table = pd.read_csv(ci['peak_area_ci_path'])
            self.assertFalse(any('1022' in c or '1018' in c or 'Threshold' in c for c in ci_table.columns))


if __name__ == '__main__':
    unittest.main()
