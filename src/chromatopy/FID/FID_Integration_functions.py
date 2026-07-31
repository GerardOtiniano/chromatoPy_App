import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy.integrate import simpson
from scipy.optimize import curve_fit
import math
from scipy.integrate import simpson
from scipy import sparse
from scipy.sparse.linalg import spsolve
from pybaselines import Baseline
import warnings
import pandas as pd
from tqdm import tqdm
from scipy.special import erf
import matplotlib.pyplot as plt

# Functions
def baseline( x, y, deg=5, max_it=1000, tol=1e-4):
    original_y = y.copy()
    order = deg + 1
    coeffs = np.ones(order)
    cond = math.pow(abs(y).max(), 1.0 / order)
    x = np.linspace(0.0, cond, y.size)  # Ensure this generates the expected range
    base = y.copy()
    vander = np.vander(x, order)  # Could potentially generate huge matrix if misconfigured
    vander_pinv = np.linalg.pinv(vander)
    for _ in range(max_it):
        coeffs_new = np.dot(vander_pinv, y)
        if np.linalg.norm(coeffs_new - coeffs) / np.linalg.norm(coeffs) < tol:
            break
        coeffs = coeffs_new
        base = np.dot(vander, coeffs)
        y = np.minimum(y, base)

    # Calculate maximum peak amplitude (3 x baseline amplitude)
    baseline_fitter = Baseline(x)
    fit, params_mask = baseline_fitter.std_distribution(y, 45)#, smooth_half_window=10)
    mask = params_mask['mask'] #  Mask for regions of signal without peaks
    min_peak_amp = (np.std(y[mask]))*2*3 # 2 sigma times 3
    return base, min_peak_amp # return base


def asls_baseline(y, lam=1e6, p=0.001, max_iter=50, conv_thresh=1e-6, return_info=True):
    """Asymmetric Least Squares baseline matching the HPLC integration path."""
    y = np.asarray(y, dtype=float).copy()
    n = y.size
    if n < 3:
        b = np.maximum(y, 0.0)
        info = {'iterations': 0, 'converged': True, 'last_delta': 0.0, 'weights': np.ones_like(y)}
        return (b, info) if return_info else b

    nan_mask = ~np.isfinite(y)
    if nan_mask.any():
        xi = np.arange(n)
        finite_mask = ~nan_mask
        if finite_mask.any():
            y[nan_mask] = np.interp(xi[nan_mask], xi[finite_mask], y[finite_mask])
        else:
            y[:] = 0.0

    diagonals = [np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)]
    offsets = [0, 1, 2]
    d_matrix = sparse.diags(diagonals, offsets, shape=(n - 2, n), format='csc')
    penalty = (d_matrix.T @ d_matrix).tocsc()

    weights = np.ones(n)
    baseline_values = y.copy()
    delta = 0.0
    for iteration in range(1, max_iter + 1):
        weight_matrix = sparse.diags(weights, 0, shape=(n, n), format='csc')
        lhs = weight_matrix + lam * penalty
        rhs = weights * y
        next_baseline = spsolve(lhs, rhs)

        residual = y - next_baseline
        weights = p * (residual > 0.0) + (1.0 - p) * (residual <= 0.0)
        weights = np.clip(weights, 1e-6, 1.0)

        denominator = np.linalg.norm(baseline_values) + 1e-12
        delta = np.linalg.norm(next_baseline - baseline_values) / denominator
        baseline_values = next_baseline
        if delta < conv_thresh:
            info = {'iterations': iteration, 'converged': True, 'last_delta': float(delta), 'weights': weights}
            return (np.maximum(baseline_values, 0.0), info) if return_info else np.maximum(baseline_values, 0.0)

    info = {'iterations': max_iter, 'converged': False, 'last_delta': float(delta), 'weights': weights}
    return (np.maximum(baseline_values, 0.0), info) if return_info else np.maximum(baseline_values, 0.0)


def hplc_style_baseline(x, y):
    """Return the same ASLS baseline and peak threshold used by HPLC integration."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx = np.median(np.diff(x)) if len(x) > 1 else 1.0
    span = (x.max() - x.min()) if len(x) else 1.0

    lam = 1e6 * max(1.0, (span / max(dx, 1e-6)) / 200.0)
    p = 0.01

    baseline_values, _ = asls_baseline(y, lam=lam, p=p, max_iter=50, conv_thresh=1e-6, return_info=True)
    baseline_values = np.maximum(baseline_values, 0.0)
    corrected = np.clip(y - baseline_values, 0, None)

    diff = np.diff(corrected)
    mad = 1.4826 * np.median(np.abs(diff - np.median(diff))) if diff.size else 0.0
    sigma = (mad / np.sqrt(2.0)) if (mad > 0 and np.isfinite(mad)) else 0.0
    dynamic_range = np.nanpercentile(y, 99) - np.nanpercentile(y, 1)
    absolute_floor = 0.005 * dynamic_range
    relative_floor = 0.02 * np.nanmedian(baseline_values) if np.isfinite(np.nanmedian(baseline_values)) else 0.0
    min_peak_amp = max(5.0 * sigma, absolute_floor, relative_floor)
    return baseline_values, float(min_peak_amp * 3)

def find_valleys(y, peaks, peak_oi=None):
    valleys = []
    if peak_oi == None:
        for i in range(1, len(peaks)):
            valley_point = np.argmin(y[peaks[i - 1] : peaks[i]]) + peaks[i - 1]
            valleys.append(valley_point)
    else:
        poi = np.where(peaks == peak_oi)[0][0]
        valleys.append(np.argmin(y[peaks[poi - 1] : peaks[poi]]) + peaks[poi - 1])
        valleys.append(np.argmin(y[peaks[poi] : peaks[poi + 1]]) + peaks[poi])
    return valleys

# def smoother(y, param_0, param_1, mode = "interp"):# "constant"):
#     return savgol_filter(y, param_0, param_1, mode=mode)
def smoother(y, window_length, polyorder):
    from scipy.signal import savgol_filter

    sample_count = len(y)
    if sample_count < 3:
        return y  # don't try to smooth tiny series
    # Ensure odd length window for Savitky-Golay
    largest_valid_window = sample_count if sample_count % 2 else sample_count - 1
    window_length = min(max(int(window_length), 1), largest_valid_window)
    if window_length % 2 == 0:
        window_length -= 1
    if window_length < 3:
        return y
    polyorder = min(max(int(polyorder), 0), window_length - 1)
    return savgol_filter(y, window_length=window_length, polyorder=polyorder)

def find_peak_neighborhood_boundaries(x, y_smooth, peaks, valleys, peak_idx,
                                      max_peaks, peak_properties, gi,
                                      smoothing_params, pk_sns):
    """Find the same transitive-overlap neighborhood used by HPLC."""
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y_smooth, dtype=float)
    peaks_arr = np.asarray(peaks, dtype=int)
    valleys_arr = np.asarray(valleys, dtype=int)
    peak_idx = int(np.clip(peak_idx, 0, max(len(x_arr) - 1, 0)))

    if x_arr.size == 0 or y_arr.size != x_arr.size:
        return None, None, []
    if peaks_arr.size == 0:
        dx = np.median(np.diff(x_arr)) if x_arr.size > 1 else 0.01
        center = float(x_arr[peak_idx])
        half_width = max(5.0 * dx, 1e-3)
        return center - half_width, center + half_width, []
    if peak_idx not in peaks_arr:
        peaks_arr = np.append(peaks_arr, peak_idx)

    peaks_arr = np.unique(peaks_arr)
    peaks_arr = peaks_arr[(peaks_arr >= 0) & (peaks_arr < len(x_arr))]
    peaks_arr = peaks_arr[np.argsort(x_arr[peaks_arr])]
    max_peaks = max(int(max_peaks or 1), 3)
    left_bases = np.asarray(peak_properties.get("left_bases", []), dtype=int)
    right_bases = np.asarray(peak_properties.get("right_bases", []), dtype=int)
    extended_boundaries = {}

    def intervals_overlap(first, second):
        first_left, first_right = sorted(first)
        second_left, second_right = sorted(second)
        return first_left <= second_right and second_left <= first_right

    def fallback_interval(peak, dx):
        center = float(x_arr[peak])
        half_width = max(5.0 * dx, 1e-3)
        return center - half_width, center + half_width

    def fit_extend_peak(peak):
        peak = int(peak)
        if peak in extended_boundaries:
            return extended_boundaries[peak]

        dx_med = np.median(np.diff(x_arr)) if len(x_arr) > 1 else 0.01
        if not np.isfinite(dx_med) or dx_med <= 0:
            dx_med = 0.01

        positions = np.flatnonzero(peaks_arr == peak)
        position = int(positions[0]) if positions.size else -1
        if (0 <= position < left_bases.size and
                position < right_bases.size):
            left = int(left_bases[position])
            right = int(right_bases[position])
        else:
            left = max(peak - 5, 0)
            right = min(peak + 5, len(x_arr) - 1)

        left_valley = max(
            (int(valley) for valley in valleys_arr if valley < peak),
            default=left,
        )
        right_valley = min(
            (int(valley) for valley in valleys_arr if valley > peak),
            default=right,
        )
        left = max(left, left_valley, 0)
        right = min(right, right_valley, len(x_arr) - 1)
        if right <= left:
            extended_boundaries[peak] = fallback_interval(peak, dx_med)
            return extended_boundaries[peak]

        sigma_estimate = max((x_arr[right] - x_arr[left]) / 20.0,
                             3.0 * dx_med)
        max_span = int(np.ceil(6.0 * sigma_estimate / dx_med))
        left = max(peak - max_span, left, 0)
        right = min(peak + max_span, right, len(x_arr) - 1)
        if right - left > 1000:
            left = max(peak - 200, 0)
            right = min(peak + 200, len(x_arr) - 1)

        x_window = x_arr[left:right + 1]
        y_window = y_arr[left:right + 1]
        if x_window.size < 3 or np.all(y_window == 0):
            extended_boundaries[peak] = fallback_interval(peak, dx_med)
            return extended_boundaries[peak]

        local_peak = int(np.argmin(np.abs(x_window - x_arr[peak])))
        try:
            heights, means, widths = estimate_initial_gaussian_params(
                pd.Series(x_window), pd.Series(y_window), local_peak)
            height = max(float(heights[0]), 0.0)
            mean = float(means[0])
            width = float(widths[0])
        except (IndexError, KeyError, TypeError, ValueError):
            height = max(float(y_arr[peak]), 0.0)
            mean = float(x_arr[peak])
            width = (x_window[-1] - x_window[0]) / 6.0

        width_max = max((x_window[-1] - x_window[0]) / 2.0, dx_med)
        width = float(np.clip(width, dx_med, width_max))
        bounds = (
            [0.0, float(x_window[0]), dx_med],
            [np.inf, float(x_window[-1]), width_max],
        )
        popt = None
        for maxfev in (gi, gi * 5):
            try:
                popt, _ = curve_fit(
                    individual_gaussian, x_window, y_window,
                    p0=[height, mean, width], bounds=bounds,
                    method="trf", maxfev=maxfev)
                break
            except (RuntimeError, ValueError):
                continue
        if popt is None:
            extended_boundaries[peak] = fallback_interval(peak, dx_med)
            return extended_boundaries[peak]

        amplitude, center, width = map(float, popt)
        x_min, x_max = calculate_gaus_extension_limits(
            center, max(width, 1e-6), factor=1)
        extended_x, extended_y = extrapolate_gaussian(
            x_window, amplitude, center, width, None,
            x_min, x_max, step=0.01)
        if extended_x.size < 3:
            extended_boundaries[peak] = fallback_interval(peak, dx_med)
            return extended_boundaries[peak]

        extended_peak = int(np.argmin(np.abs(extended_x - center)))
        left_idx, right_idx = calculate_hplc_boundaries(
            extended_x, extended_y, extended_peak, smoothing_params)
        extended_boundaries[peak] = (
            float(extended_x[left_idx]), float(extended_x[right_idx]))
        return extended_boundaries[peak]

    peak_position = int(np.flatnonzero(peaks_arr == peak_idx)[0])
    included = {peak_idx}
    left_position = right_position = peak_position

    while left_position > 0 and len(included) < max_peaks:
        current_peak = int(peaks_arr[left_position])
        left_peak = int(peaks_arr[left_position - 1])
        if not intervals_overlap(
                fit_extend_peak(current_peak), fit_extend_peak(left_peak)):
            break
        included.add(left_peak)
        left_position -= 1

    while right_position < len(peaks_arr) - 1 and len(included) < max_peaks:
        current_peak = int(peaks_arr[right_position])
        right_peak = int(peaks_arr[right_position + 1])
        if not intervals_overlap(
                fit_extend_peak(current_peak), fit_extend_peak(right_peak)):
            break
        included.add(right_peak)
        right_position += 1

    boundaries = [fit_extend_peak(peak) for peak in included]
    neighborhood_left = min(boundary[0] for boundary in boundaries)
    neighborhood_right = max(boundary[1] for boundary in boundaries)
    overlapping_peaks = sorted(peak for peak in included if peak != peak_idx)
    return neighborhood_left, neighborhood_right, overlapping_peaks


# Gaussian fitting
# def calculate_gaus_extension_limits(cen, wid, decay, factor=2, max_tail_sigma=2):#5):
#     sigma_effective = wid * factor  # Adjust factor for tail thinness
#     if decay <= 0:
#         tail = sigma_effective * max_tail_sigma
#     else:
#         tail = min(1/decay, sigma_effective * max_tail_sigma)
#     return cen - sigma_effective-tail, cen+sigma_effective+tail
def calculate_gaus_extension_limits(cen, wid, factor=2, max_tail_sigma=3):
    sigma_effective = wid * factor
    tail = sigma_effective * max_tail_sigma
    return cen - sigma_effective - tail, cen + sigma_effective + tail

def extrapolate_gaussian(x, amp, cen, wid, skew=None, x_min=None, x_max=None, step=0.001):
    if x_min is None: x_min = cen - 3 * wid
    if x_max is None: x_max = cen + 3 * wid
    extended_x = np.arange(x_min, x_max, step)
    if skew is None:
        extended_y = individual_gaussian(extended_x, amp, cen, wid)
    else:
        extended_y = skewed_gaussian(extended_x, amp, cen, wid, skew)
    return extended_x, extended_y

# def extrapolate_gaussian_decay(amp, cen, wid, dec, x_min=None, x_max=None, step=1e-4):
#     if x_min is None:
#         x_min = cen - 3 * wid
#     if x_max is None:
#         x_max = cen + 3 * wid
#     xs = np.arange(x_min, x_max, step)
#     ys = gaussian_decay(xs, amp, cen, wid, dec)
#     return xs, ys

def calculate_boundaries(x, y, ind_peak, smoothing_params, pk_sns):
    smooth_y = smoother(y, smoothing_params[0], smoothing_params[1])
    velocity, X1 = forward_derivative(x, smooth_y)
    velocity /= np.max(np.abs(velocity))
    if smoothing_params[0] > len(velocity):
        smoother_val = len(velocity)-1
    else: smoother_val = smoothing_params[0]
    smooth_velo = smoother(velocity, smoother_val, smoothing_params[1])
    dt = int(np.ceil(0.025 / np.mean(np.diff(x))))
    A = np.where(smooth_velo[: ind_peak - 3 * dt] < pk_sns)[0]  # 0.05)[0]
    B = np.where(smooth_velo[ind_peak + 3 * dt :] > -pk_sns)[0]  # -0.05)[0]
    if A.size > 0:
        A = A[-1] + 1
    else:
        A = 1
    if B.size > 0:
        B = B[0] + ind_peak + 3 * dt - 1
    else:
        B = len(x) - 1
    return A, B


def calculate_hplc_boundaries(x, y, ind_peak, smoothing_params,
                              tolerance=0.02, w_factor=3.0):
    """HPLC's local derivative boundary finder, used by MultiGaussian fits."""
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if x_arr.size < 3 or y_arr.size != x_arr.size:
        return 0, max(len(x_arr) - 1, 0)

    y_smooth = np.asarray(
        smoother(y_arr, smoothing_params[0], smoothing_params[1]),
        dtype=float,
    )
    velocity, _ = forward_derivative(x_arr, y_smooth)
    velocity = np.asarray(velocity, dtype=float)
    if velocity.size == x_arr.size - 1 and velocity.size:
        velocity = np.r_[velocity, velocity[-1]]
    if velocity.size != x_arr.size:
        velocity = np.resize(velocity, x_arr.size)
    velocity[~np.isfinite(velocity)] = 0.0
    velocity_smooth = np.asarray(
        smoother(velocity, smoothing_params[0], smoothing_params[1]),
        dtype=float,
    )
    absolute_velocity = np.abs(velocity_smooth)

    peak = int(np.clip(ind_peak, 0, len(y_smooth) - 1))
    half_height = 0.5 * y_smooth[peak]
    left = peak
    while left > 0 and y_smooth[left] > half_height:
        left -= 1
    right = peak
    while right < len(y_smooth) - 1 and y_smooth[right] > half_height:
        right += 1

    window = max(5, int((right - left) * w_factor / 2.0))
    local_left = max(0, peak - window)
    local_right = min(len(y_smooth) - 1, peak + window)
    if peak > local_left:
        left_shoulder = local_left + int(
            np.argmax(absolute_velocity[local_left:peak]))
    else:
        left_shoulder = peak
    if local_right > peak + 1:
        right_shoulder = peak + 1 + int(
            np.argmax(absolute_velocity[peak + 1:local_right + 1]))
    else:
        right_shoulder = peak

    epsilon = 1e-12
    left_threshold = max(
        tolerance * abs(velocity_smooth[left_shoulder]), epsilon)
    right_threshold = max(
        tolerance * abs(velocity_smooth[right_shoulder]), epsilon)
    left_boundary = left_shoulder
    while (left_boundary > local_left and
           absolute_velocity[left_boundary] > left_threshold):
        left_boundary -= 1
    right_boundary = right_shoulder
    while (right_boundary < local_right and
           absolute_velocity[right_boundary] > right_threshold):
        right_boundary += 1

    if left_boundary >= right_boundary:
        left_boundary = max(local_left, peak - 3)
        right_boundary = min(local_right, peak + 3)
    return int(left_boundary), int(right_boundary)


def calculate_asymmetric_boundaries(y, relative_height=0.001):
    """Return model-aware boundaries for a fitted asymmetric peak.

    The skew-normal fit already contains the peak asymmetry, so its left and
    right limits should be found independently.  Retain the fitted curve until
    both tails fall to the same fraction of the fitted maximum instead of
    imposing a symmetric search window around the apex.
    """
    y_arr = np.asarray(y, dtype=float)
    if y_arr.size == 0:
        return 0, 0
    if not 0 < relative_height < 1:
        raise ValueError("relative_height must be between 0 and 1")

    finite_y = np.where(np.isfinite(y_arr), np.maximum(y_arr, 0.0), 0.0)
    peak_height = float(np.max(finite_y))
    if peak_height <= 0:
        return 0, len(finite_y) - 1

    retained = np.flatnonzero(finite_y >= relative_height * peak_height)
    if retained.size == 0:
        peak = int(np.argmax(finite_y))
        return peak, peak
    return int(retained[0]), int(retained[-1])


def calculate_boundaries_acceleration(x, y, ind_peak, smoothing_params, pk_sns):
    smooth_y = smoother(y, smoothing_params[0], smoothing_params[1])
    velocity, _ = forward_derivative(x, smooth_y)
    acceleration, _ = forward_derivative(x[:-1], velocity)
    acceleration /= np.max(np.abs(acceleration))
    smoother_val = min(smoothing_params[0], len(acceleration) - 1)
    smooth_accel = smoother(acceleration, smoother_val, smoothing_params[1])
    left_zone = smooth_accel[:ind_peak]
    right_zone = smooth_accel[ind_peak:]
    if len(left_zone) > 0:
        A = np.argmax(left_zone)
    else:
        A = 1
    if len(right_zone) > 0:
        B = np.argmax(right_zone) + ind_peak
    else:
        B = len(x) - 1
    return A, B
# def calculate_boundaries_acceleration(x, y, ind_peak, smoothing_params, pk_sns):
#     smooth_y = smoother(y, smoothing_params[0], smoothing_params[1])

#     # 1) Derivatives
#     vel, _ = forward_derivative(x, smooth_y)               # len = n-1
#     if vel.size == 0:
#         return 0, max(len(x) - 1, 0)

#     acc, _ = forward_derivative(x[:-1], vel)               # len = n-2
#     if acc.size == 0:
#         return 0, max(len(x) - 1, 0)

#     # 2) Normalize safely (avoid div-by-zero)
#     amax = np.max(np.abs(acc))
#     if amax > 0:
#         acc = acc / amax

#     # 3) Smoothing window must be valid (odd & >= polyorder+2)
#     win = min(smoothing_params[0], len(acc) - 1)
#     if win % 2 == 0:
#         win -= 1
#     win = max(win, smoothing_params[1] + 2 + (smoothing_params[1] % 2))
#     win = min(win, max(len(acc) - 1, 1))
#     if win < 3:  # too short to smooth; skip smoothing
#         smooth_acc = acc
#     else:
#         smooth_acc = smoother(acc, win, smoothing_params[1])

#     # 4) Split around the peak index (ind_peak refers to x/y index)
#     left_zone  = smooth_acc[:max(ind_peak, 0)]
#     right_zone = smooth_acc[max(ind_peak, 0):]

#     # 5) Pick extremum on each side; clamp to valid [0, len(x)-1]
#     A = int(np.argmax(left_zone)) if left_zone.size > 0 else 0
#     B = int(np.argmax(right_zone)) + max(ind_peak, 0) if right_zone.size > 0 else len(x) - 1

#     # 6) Ensure A < B and avoid A-1 underflow in downstream slicing
#     A = max(A, 1)                  # so slice [A-1: ...] doesn’t go negative
#     B = min(B, len(x) - 1)

#     if B <= A:
#         # conservative fallback: symmetric window around the peak
#         pad = max(int(0.02 / max(np.mean(np.diff(x)), 1e-9)), 10)  # ~window in samples
#         A = max(ind_peak - pad, 1)
#         B = min(ind_peak + pad, len(x) - 1)

#     return A, B

def fit_gaussians(x_full, y_full, ind_peak, peaks, smoothing_params, pk_sns,
                  gi, mode="both", valleys=None):
    if mode not in {"single", "multi", "both", "asymmetric", "asymmetric_or_multi"}:
        raise ValueError("mode must be 'single', 'multi', 'both', 'asymmetric', or 'asymmetric_or_multi'")
    # figy = plt.figure()
    results = []

    # --- MULTI-GAUSSIAN ---
    if mode in {"multi", "both", "asymmetric_or_multi"}:
        result = _fit_multi_gaussian(
            x_full, y_full, ind_peak, peaks, smoothing_params, pk_sns, gi,
            valleys=valleys)
        if result is not None:
            (best_x, best_fit_y, best_fit_params, best_fit_params_error,
             best_error, best_idx_interest, result_name,
             result_multi_flag) = result
            results.append({
                "name": result_name,
                "x": best_x,
                "y": best_fit_y,
                "params": best_fit_params,
                "pcov": best_fit_params_error,
                "error": best_error,
                "idx_interest": best_idx_interest,
                "multi_flag": result_multi_flag})

    # --- SINGLE-GAUSSIAN ---
    if mode in {"single", "both"}:
        result = _fit_single_gaussian(x_full, y_full, ind_peak, smoothing_params, pk_sns, gi, current_best_error=float("inf"))
        if result is not None:
            best_x, best_fit_y, best_fit_params, best_fit_params_error, best_error = result
            results.append({
                "name": "single",
                "x": best_x,
                "y": best_fit_y,
                "params": best_fit_params,
                "pcov": best_fit_params_error,
                "error": best_error,
                "multi_flag": False,
                "idx_interest": None})

    # --- ASYMMETRIC MODEL ---
    if mode in {"asymmetric", "asymmetric_or_multi", "both"}:
        result = _fit_asymmetric_gaussian(x_full, y_full, ind_peak, smoothing_params, pk_sns, gi, current_best_error=float("inf"))
        if result is not None:
            best_x, best_fit_y, best_fit_params, best_fit_params_error, best_error = result
            results.append({
                "name": "asymmetric",
                "x": best_x,
                "y": best_fit_y,
                "params": best_fit_params,
                "pcov": best_fit_params_error,
                "error": best_error,
                "multi_flag": False,
                "idx_interest": None})
    if not results:
        raise RuntimeError(f"No valid fit found for peak at index {ind_peak}")

    best_result = min(results, key=lambda r: r["error"])
    # --- Process best fit output ---
    best_x = best_result["x"]
    best_fit_y = best_result["y"]
    best_fit_params = best_result["params"]
    best_fit_params_error = best_result["pcov"]
    best_idx_interest = best_result.get("idx_interest", None)
    multi_gauss_flag = best_result["multi_flag"]
    model_used = best_result["name"]
    integration_x = None
    integration_y = None
    # Extend the selected component from its fitted width (same as HPLC)
    if multi_gauss_flag:
        amp, cen, wid = best_fit_params[best_idx_interest * 3: best_idx_interest * 3 + 3]
        x_min, x_max = calculate_gaus_extension_limits(cen, wid)
        best_x, best_fit_y = extrapolate_gaussian(
            best_x, amp, cen, wid, None, x_min, x_max, step=0.001)
        new_ind_peak = (np.abs(best_x - x_full[ind_peak])).argmin()
        left_boundary, right_boundary = calculate_hplc_boundaries(
            best_x, best_fit_y, new_ind_peak, smoothing_params)
        segment_x, segment_y = _clamped_segment(best_x, best_fit_y, left_boundary, right_boundary)
        if segment_x is not None:
            best_x, best_fit_y = segment_x, segment_y
    else:
        amp, cen, wid = best_fit_params[:3]
        tail_factor = 2
        if model_used == "asymmetric":
            alpha = best_fit_params[3]
            x_min, x_max = calculate_gaus_extension_limits(cen, wid, factor=tail_factor)
            best_x, best_fit_y = extrapolate_gaussian(
                best_x, amp, cen, wid, alpha, x_min, x_max, step=0.001)
            # Integrate the complete fitted curve.  The display boundary below
            # is intentionally independent so plotting cannot truncate area.
            integration_x = best_x
            integration_y = best_fit_y
        else:
            x_min, x_max = calculate_gaus_extension_limits(cen, wid, factor=tail_factor)
            best_x, best_fit_y = extrapolate_gaussian(
                best_x, amp, cen, wid, None, x_min, x_max, step=0.001)
        if model_used == "asymmetric":
            left_boundary, right_boundary = calculate_asymmetric_boundaries(
                best_fit_y, relative_height=0.001)
        elif model_used == "single":
            new_ind_peak = (np.abs(best_x - x_full[ind_peak])).argmin()
            left_boundary, right_boundary = calculate_hplc_boundaries(
                best_x, best_fit_y, new_ind_peak, smoothing_params)
        else:
            new_ind_peak = (np.abs(best_x - x_full[ind_peak])).argmin()
            left_boundary, right_boundary = calculate_boundaries_acceleration(
                best_x, best_fit_y, new_ind_peak, smoothing_params, pk_sns)
        segment_x, segment_y = _clamped_segment(best_x, best_fit_y, left_boundary, right_boundary)
        if segment_x is not None:
            best_x, best_fit_y = segment_x, segment_y

    area_x = integration_x if integration_x is not None else best_x
    area_y = integration_y if integration_y is not None else best_fit_y
    area_smooth = float(simpson(y=np.maximum(area_y, 0.0), x=area_x)) if len(area_x) >= 2 else 0.0
    return best_x, best_fit_y, area_smooth, _model_result_for_storage(best_result)


def _model_result_for_storage(best_result):
    """Return selected and full fit metadata for post-calculated uncertainty."""
    all_params = np.asarray(best_result["params"], dtype=float)
    all_covariance = np.asarray(best_result["pcov"], dtype=float)
    model_name = best_result["name"]

    if best_result["multi_flag"]:
        selected_index = int(best_result["idx_interest"])
        start, end = 3 * selected_index, 3 * selected_index + 3
        parameter_names = ["Amplitude", "Center", "Width"]
        selected_params = all_params[start:end]
        selected_covariance = all_covariance[start:end, start:end]
    else:
        selected_index = 0
        parameter_names = ["Amplitude", "Center", "Width"]
        if model_name == "asymmetric":
            parameter_names.append("Alpha")
        n_params = len(parameter_names)
        selected_params = all_params[:n_params]
        selected_covariance = all_covariance[:n_params, :n_params]

    uncertainties = np.sqrt(np.maximum(np.diag(selected_covariance), 0.0))
    model_types = {
        "multi": "Multi-Gaussian",
        "single": "Single-Gaussian",
        "asymmetric": "Asymmetric-Gaussian",
    }
    return {
        "name": model_name,
        "model_type": model_types[model_name],
        "parameter_names": parameter_names,
        "parameters": selected_params,
        "parameter_uncertainties": uncertainties,
        "covariance": selected_covariance,
        "all_parameters": all_params,
        "all_covariance": all_covariance,
        "selected_component_index": selected_index,
        "error": float(best_result["error"]),
    }


def build_processed_peak_result(area, model_result, retention_time, fit_x, fit_y):
    """Format a persisted FID peak like HPLC, plus covariance metadata."""
    params = {
        name: float(value)
        for name, value in zip(model_result["parameter_names"], model_result["parameters"])
    }
    params.update({
        f"{name} Unc": float(value)
        for name, value in zip(
            model_result["parameter_names"], model_result["parameter_uncertainties"])
    })
    params.update({
        "Covariance": model_result["covariance"],
        "All Parameters": model_result["all_parameters"],
        "All Parameter Covariance": model_result["all_covariance"],
        "Selected Component Index": model_result["selected_component_index"],
    })
    return {
        "Peak Area - best fit": float(area),
        "Retention Time": float(retention_time),
        "Model Type": model_result["model_type"],
        "Model Parameters": params,
        "Fit Error": model_result["error"],
        "Fit": {
            "x": np.asarray(fit_x, dtype=float),
            "y": np.asarray(fit_y, dtype=float),
        },
    }


def _sigma_from_curvature(x_values, y_values, local_peak, smoothing_params,
                          epsilon=1e-12):
    """Estimate HPLC's local Gaussian width from apex curvature."""
    x_arr = np.asarray(x_values, dtype=float)
    y_arr = np.asarray(y_values, dtype=float)
    if local_peak < 1 or local_peak > len(x_arr) - 2:
        return None
    y_smooth = np.asarray(
        smoother(y_arr, smoothing_params[0], smoothing_params[1]),
        dtype=float,
    )
    dx_left = x_arr[local_peak] - x_arr[local_peak - 1]
    dx_right = x_arr[local_peak + 1] - x_arr[local_peak]
    second_derivative = 2.0 * (
        (y_smooth[local_peak + 1] - y_smooth[local_peak]) /
        (dx_right + epsilon)
        - (y_smooth[local_peak] - y_smooth[local_peak - 1]) /
        (dx_left + epsilon)
    ) / (dx_left + dx_right + epsilon)
    amplitude = max(float(y_smooth[local_peak]), epsilon)
    if not np.isfinite(second_derivative) or second_derivative >= -epsilon:
        return None
    return float(np.sqrt(amplitude / (-second_derivative + epsilon)))


def _width_seed_and_bounds(x_segment, y_segment, center, left_center,
                           right_center, width_min, width_max,
                           smoothing_params, gap_factor=0.6):
    """Apply HPLC's curvature seed and nearest-neighbor width cap."""
    x_arr = np.asarray(x_segment, dtype=float)
    y_arr = np.asarray(y_segment, dtype=float)
    local_peak = int(np.argmin(np.abs(x_arr - center)))
    curvature_width = _sigma_from_curvature(
        x_arr, y_arr, local_peak, smoothing_params)
    left_gap = abs(center - left_center) if left_center is not None else np.inf
    right_gap = abs(center - right_center) if right_center is not None else np.inf
    nearest_gap = min(left_gap, right_gap)
    gap_cap = gap_factor * nearest_gap if np.isfinite(nearest_gap) else np.inf
    upper = min(gap_cap, width_max) if np.isfinite(gap_cap) else width_max
    upper = max(float(upper), float(width_min))
    if curvature_width is None or not np.isfinite(curvature_width):
        initial = upper
    else:
        initial = min(float(curvature_width), upper)
    initial = float(np.clip(initial, width_min, upper))
    return initial, float(width_min), upper


def _closer_boundary(x_arr, derivative_boundary, valley_boundary, apex):
    if valley_boundary is None:
        return int(derivative_boundary)
    derivative_distance = abs(float(x_arr[derivative_boundary]) - apex)
    valley_distance = abs(float(x_arr[valley_boundary]) - apex)
    return (int(derivative_boundary) if derivative_distance <= valley_distance
            else int(valley_boundary))


def _fit_hplc_single_fallback(x_arr, y_arr, ind_peak, smoothing_params,
                              valleys, gi):
    """Fit the tightly bounded single-Gaussian fallback from HPLC."""
    left, right = calculate_hplc_boundaries(
        x_arr, y_arr, ind_peak, smoothing_params)
    apex = float(x_arr[ind_peak])
    valleys_arr = np.asarray(valleys if valleys is not None else [], dtype=int)
    valleys_arr = valleys_arr[(valleys_arr >= 0) & (valleys_arr < len(x_arr))]
    left_valleys = valleys_arr[valleys_arr < ind_peak]
    right_valleys = valleys_arr[valleys_arr > ind_peak]
    left_valley = int(left_valleys[-1]) if left_valleys.size else None
    right_valley = int(right_valleys[0]) if right_valleys.size else None
    left = _closer_boundary(x_arr, left, left_valley, apex)
    right = _closer_boundary(x_arr, right, right_valley, apex)
    if not (0 <= left < ind_peak < right < len(x_arr)):
        left, right = calculate_hplc_boundaries(
            x_arr, y_arr, ind_peak, smoothing_params)
    if right - left < 2:
        return None

    x_segment = x_arr[left:right + 1]
    y_segment = y_arr[left:right + 1]
    local_peak = int(np.argmin(np.abs(x_segment - apex)))
    heights, _, widths = estimate_initial_gaussian_params(
        pd.Series(x_segment), pd.Series(y_segment), local_peak)
    amplitude = max(float(heights[0]), 0.0)
    width = float(widths[0])
    differences = np.diff(x_segment)
    finite_differences = differences[
        np.isfinite(differences) & (differences > 0)]
    dx = float(np.median(finite_differences)) if finite_differences.size else 1.0
    width_min = max(1.5 * dx, 1e-3)
    width_max = max((x_segment[-1] - x_segment[0]) / 3.0,
                    2.0 * width_min)
    width = float(np.clip(width, width_min, width_max))
    amplitude_max = max(1.0 + float(y_arr[ind_peak]), amplitude * 3.0)
    bounds = (
        [0.0, apex - 0.01, width_min],
        [amplitude_max, apex + 0.01, width_max],
    )
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(
                individual_gaussian, x_segment, y_segment,
                p0=[amplitude, apex, width], bounds=bounds,
                method="trf", maxfev=gi)
    except (RuntimeError, ValueError):
        return None
    fitted = individual_gaussian(x_segment, *popt)
    error = float(np.sqrt(np.mean((fitted - y_segment) ** 2)))
    return x_segment, fitted, popt, pcov, error


def _fit_multi_gaussian(x_full, y_full, ind_peak, peaks, smoothing_params,
                        pk_sns, gi, valleys=None):
    """Run HPLC's iterative MultiGaussian fit and single-peak fallback."""
    x_arr = np.asarray(x_full, dtype=float)
    y_arr = np.asarray(y_full, dtype=float)
    ind_peak = int(ind_peak)
    current_peaks = np.unique(
        np.append(np.asarray(peaks, dtype=int), ind_peak))
    current_peaks = current_peaks[
        (current_peaks >= 0) & (current_peaks < len(x_arr))]
    current_peaks = current_peaks[np.argsort(x_arr[current_peaks])]
    if ind_peak not in current_peaks or x_arr.size < 3:
        return None

    best_fit_y = None
    best_fit_params = None
    best_fit_covariance = None
    best_x = None
    best_error = float("inf")
    best_idx_interest = None
    best_name = "multi"
    best_multi_flag = True
    full_dx = np.diff(x_arr)
    full_dx = full_dx[np.isfinite(full_dx) & (full_dx > 0)]
    dx_full = float(np.median(full_dx)) if full_dx.size else 1.0

    while len(current_peaks) > 1:
        left, _ = calculate_hplc_boundaries(
            x_arr, y_arr, int(current_peaks[0]), smoothing_params)
        _, right = calculate_hplc_boundaries(
            x_arr, y_arr, int(current_peaks[-1]), smoothing_params)
        left = max(int(left), 0)
        right = min(int(right), len(x_arr) - 1)
        if right - left < 2:
            break
        x_segment = x_arr[left:right + 1]
        y_segment = y_arr[left:right + 1]

        current_peaks = current_peaks[np.argsort(x_arr[current_peaks])]
        seed_centers = x_arr[current_peaks]
        left_neighbors = [
            seed_centers[index - 1] if index > 0 else None
            for index in range(len(seed_centers))]
        right_neighbors = [
            seed_centers[index + 1] if index + 1 < len(seed_centers) else None
            for index in range(len(seed_centers))]

        segment_dx = np.diff(x_segment)
        segment_dx = segment_dx[
            np.isfinite(segment_dx) & (segment_dx > 0)]
        dx = float(np.median(segment_dx)) if segment_dx.size else dx_full
        width_min = max(1.5 * dx, 1e-3)
        width_max = max((x_segment[-1] - x_segment[0]) / 3.0,
                        2.0 * width_min)
        initial = []
        lower = []
        upper = []
        for index, peak in enumerate(current_peaks):
            amplitude = max(float(y_arr[peak]), 0.0)
            center = float(x_arr[peak])
            width, width_lower, width_upper = _width_seed_and_bounds(
                x_segment, y_segment, center,
                left_neighbors[index], right_neighbors[index],
                width_min, width_max, smoothing_params,
                gap_factor=0.6)
            amplitude_upper = max(
                1.0 + float(y_arr[peak]), amplitude * 3.0)
            initial.extend([amplitude, center, width])
            lower.extend([0.0, center - 0.15, width_lower])
            upper.extend([amplitude_upper, center + 0.15, width_upper])

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, pcov = curve_fit(
                    multigaussian, x_segment, y_segment, p0=initial,
                    bounds=(lower, upper), method="trf", maxfev=gi)
            fitted = multigaussian(x_segment, *popt)
            fitted_centers = np.asarray(popt[1::3], dtype=float)
            selected_component = int(np.argmin(
                np.abs(fitted_centers - x_arr[ind_peak])))
            error = float(np.sqrt(np.mean((fitted - y_segment) ** 2)))
            if error < best_error:
                best_fit_y = fitted
                best_fit_params = popt
                best_fit_covariance = pcov
                best_x = x_segment
                best_error = error
                best_idx_interest = selected_component
        except (RuntimeError, ValueError):
            pass

        distances = np.abs(x_arr[current_peaks] - x_arr[ind_peak])
        current_peaks = np.delete(current_peaks, int(np.argmax(distances)))

    single_result = _fit_hplc_single_fallback(
        x_arr, y_arr, ind_peak, smoothing_params, valleys, gi)
    if single_result is not None:
        single_x, single_y, single_params, single_covariance, single_error = single_result
        if best_fit_params is None or single_error < best_error / 1.02:
            best_x = single_x
            best_fit_y = single_y
            best_fit_params = single_params
            best_fit_covariance = single_covariance
            best_error = single_error
            best_idx_interest = 0
            best_name = "single"
            best_multi_flag = False

    if best_fit_params is None:
        return None
    return (best_x, best_fit_y, best_fit_params, best_fit_covariance,
            best_error, best_idx_interest, best_name, best_multi_flag)

def _fit_single_gaussian(x_full, y_full, ind_peak, smoothing_params, pk_sns, gi, current_best_error):
    left, right = calculate_boundaries(x_full, y_full, ind_peak, smoothing_params, pk_sns)
    x = x_full[left:right + 1]
    y = y_full[left:right + 1]
    h, c, w = estimate_initial_gaussian_params(x, y, ind_peak)
    center_idx = (np.abs(x - c[0])).argmin()
    # decay_init = estimate_initial_decay(x, y, center_idx)
    p0 = [h[0], c[0], w[0]]#, decay_init]
    # p0 = [h[0], c[0], w[0], 0.1]
    bounds = ([0.9 * y_full[ind_peak], x_full[ind_peak] - 0.1, 0.5 * w[0]],
              [1 + y_full[ind_peak], x_full[ind_peak] + 0.1, 1.5 * w[0]])

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # popt, pcov = curve_fit(gaussian_decay, x, y, p0=p0, method="trf", bounds=bounds, maxfev=gi)
            popt, pcov = curve_fit(individual_gaussian, x, y, p0=p0, method="trf", bounds=bounds, maxfev=gi)
        # fitted_y = gaussian_decay(x, *popt)
        fitted_y = individual_gaussian(x, *popt)
        error = np.sqrt(np.mean((fitted_y - y) ** 2))
        if error < current_best_error:
            return x, fitted_y, popt, pcov, error
    except RuntimeError:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # popt, pcov = curve_fit(gaussian_decay, x, y, p0=p0, method="trf", bounds=bounds, maxfev=gi * 5)
                popt, pcov = curve_fit(individual_gaussian, x, y, p0=p0, method="trf", bounds=bounds, maxfev=gi * 5)
            # fitted_y = gaussian_decay(x, *popt)
            fitted_y = individual_gaussian(x, *popt)
            error = np.sqrt(np.mean((fitted_y - y) ** 2))
            if error < current_best_error:
                return x, fitted_y, popt, pcov, error
        except RuntimeError:
            tqdm.write("Error: Optimal parameters could not be found even after increasing iterations.")
    return None

def _fit_asymmetric_gaussian(x_full, y_full, ind_peak, smoothing_params, pk_sns, gi, current_best_error):
    left, right = calculate_boundaries(x_full, y_full, ind_peak, smoothing_params, pk_sns)
    x = x_full[left:right + 1]
    y = y_full[left:right + 1]

    h, c, w = estimate_initial_gaussian_params(x, y, ind_peak)

    # Sanitize estimates
    amp = max(h[0], 1e-5)
    cen = c[0]
    wid = max(w[0], 1e-5)
    alpha = 0.0  # Start symmetric

    p0 = [amp, cen, wid, alpha]
    bounds = (
        [1e-5, cen - 0.1, 1e-5, -10],  # lower bounds
        [10 * amp, cen + 0.1, 10 * wid, 10]  # upper bounds
    )

    try:
        popt, pcov = curve_fit(skewed_gaussian, x, y, p0=p0, method="trf", bounds=bounds, maxfev=gi)
        fitted_y = skewed_gaussian(x, *popt)
        error = np.sqrt(np.mean((fitted_y - y) ** 2))
        if error < current_best_error:
            return x, fitted_y, popt, pcov, error
    except RuntimeError:
        pass
    return None

def _clamped_segment(x_arr, y_arr, L, R):
    n = len(x_arr)
    # clamp
    L = max(L, 0)
    R = min(R, n - 1)
    if R < L:
        return None, None
    # expand by 1 safely
    L = max(L - 1, 0)
    R = min(R + 1, n - 1)
    xs = x_arr[L:R + 1]
    ys = y_arr[L:R + 1]
    if xs.size == 0 or ys.size == 0:
        return None, None
    return xs, ys
def debug_param_distribution(mu, cov, n_draw=5000):
    """
    Sample from N(mu, cov) and plot:
      1) joint scatter of (wid, decay)
      2) histogram of wid
      3) histogram of decay
    """
    mu = np.asarray(mu)
    cov = np.asarray(cov)

    # draw a big batch
    batch = np.random.multivariate_normal(mu, cov, size=n_draw)
    wid   = batch[:,2]
    decay = batch[:,3]

    # 1) Joint scatter
    plt.figure()
    plt.scatter(wid, decay, alpha=0.2)
    plt.axvline(0)
    plt.axhline(0)
    plt.xlabel("wid")
    plt.ylabel("decay")
    plt.title("Joint draw of (wid, decay)")
    plt.show()

    # 2) wid histogram
    plt.figure()
    plt.hist(wid, bins=50)
    plt.xlabel("wid")
    plt.title("Histogram of wid")
    plt.show()

    # 3) decay histogram
    plt.figure()
    plt.hist(decay, bins=50)
    plt.xlabel("decay")
    plt.title("Histogram of decay")
    plt.show()

def individual_gaussian( x, amp, cen, wid):
    return amp * np.exp(-((x - cen) ** 2) / (2 * wid**2))

def estimate_initial_gaussian_params( x, y, peak):
    # Subset peaks so that only idx positions with x bounds are considered
    heights = []
    means = []
    stddevs = []
    height = y[peak]
    mean = x[peak]
    half_max = 0.5 * height
    mask = y >= half_max
    valid_x = x[mask]
    if len(valid_x) > 1:
        fwhm = np.abs(valid_x.iloc[-1] - valid_x.iloc[0])
        stddev = fwhm / (2 * np.sqrt(2 * np.log(2)))
    else:
        stddev = (x.max() - x.min()) / 6
    heights.append(height)
    means.append(mean)
    stddevs.append(stddev)
    return heights, means, stddevs

def estimate_initial_decay(x, y, center_idx):
    left_half = y[:center_idx]
    right_half = y[center_idx:]
    left_slope = np.mean(np.gradient(left_half))
    right_slope = np.mean(np.gradient(right_half))
    asymmetry = right_slope - left_slope

    # Empirical mapping to decay (tweak this based on real data behavior)
    decay_est = np.clip(0.5 * asymmetry, 0.01, 1.5)
    return decay_est

def multigaussian( x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amp = params[i]
        cen = params[i + 1]
        wid = params[i + 2]
        y += amp * np.exp(-((x - cen) ** 2) / (2 * wid**2))
    return y

def skewed_gaussian(x, amp, cen, sigma, alpha):
    """
    Skewed Gaussian (Skew-Normal) distribution:
    - alpha = 0 gives symmetric Gaussian
    - alpha > 0 → right skew
    - alpha < 0 → left skew
    """
    z = (x - cen) / (sigma * np.sqrt(2))
    return amp * np.exp(-z**2) * (1 + erf(alpha * z))

def gaussian_decay( x, amp, cen, wid, dec):
    return amp * np.exp(-((x - cen) ** 2) / (2 * wid**2)) * np.exp(-dec * abs(x - cen))


def fit_color(model_params):
    model_name = ""
    if isinstance(model_params, dict):
        model_name = str(model_params.get("name", "")).lower()
    return "blue" if model_name == "asymmetric" else "red"


def add_peak_label(ax, label, x_fit, y_fit, x_signal, y_signal):
    if len(x_fit) == 0 or len(y_fit) == 0:
        return None
    fit_peak_idx = int(np.argmax(y_fit))
    x_peak = float(x_fit[fit_peak_idx])
    y_peak = float(y_fit[fit_peak_idx])
    y_values = np.asarray(y_signal, dtype=float)
    y_min, y_max = ax.get_ylim()
    y_range = max(y_max - y_min, np.nanmax(y_values) - np.nanmin(y_values), 1e-9)

    x_values = np.asarray(x_signal, dtype=float)
    local_mask = np.abs(x_values - x_peak) <= 0.08
    if np.any(local_mask):
        y_peak = max(y_peak, float(np.nanmax(y_values[local_mask])))

    y_text = y_peak + 0.05 * y_range
    if y_text > y_max - 0.04 * y_range:
        ax.set_ylim(y_min, y_text + 0.08 * y_range)

    return ax.text(
        x_peak,
        y_text,
        str(label),
        ha="center",
        va="bottom",
        fontsize=8,
        color="black",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.65, pad=1.5),
        zorder=3)

def forward_derivative(x, y):
    fd = np.diff(y) / np.diff(x)
    x_n = x#[:-1]
    return fd, x_n

class FIDAnalyzer:
    def __init__(self, df, window_bounds, gaus_iterations, sample_name, is_reference, max_peaks, sw, sf, pk_sns, pk_pr, max_PA, reference_peaks=None):
        self.fig, self.axs = None, None
        self.df = df
        self.window_bounds = window_bounds
        self.sample_name = sample_name
        self.is_reference = is_reference
        self.reference_peaks = reference_peaks  # ref_key
        self.fig, self.axs = None, None
        self.datasets = []
        self.peaks_indices = []
        self.integrated_peaks = {}
        self.action_stack = []
        self.no_peak_lines = {}
        self.peaks = {}  # Store all peak indices and properties for each trace
        self.axs_to_traces = {}  # Empty map for connecting traces to figure axes
        self.peak_results = {}
        self.peak_results['Sample ID'] = sample_name
        self.gi = gaus_iterations
        self.max_peaks_for_neighborhood = max_peaks
        self.peak_properties = {}
        self.smoothing_params = [sw, sf]
        self.pk_sns = pk_sns
        self.pk_pr = pk_pr
        self.t_pressed = False # Flag to track if 't' was pressed
        self.called = False
        self.max_peak_amp = max_PA

    def run(self):
        """
        Executes the peak analysis workflow.
        Returns:
            peaks (dict): Peak areas and related info.
            fig (matplotlib.figure.Figure): The figure object.
            reference_peaks (dict): Updated reference peaks.
            t_pressed (bool): Indicates if 't' was pressed to update reference peaks.
        """
        self.fig, self.axs = self.plot_data()
        self.current_ax_idx = 0  # Initialize current axis index
        if self.is_reference:
            # Reference samples handling
            self.fig.canvas.mpl_connect("button_press_event", self.on_click)
            self.fig.canvas.mpl_connect("key_press_event", self.on_key)  # Connect general key events
            plt.show(block=True)  # Blocks script until plot window is closed
            if not self.reference_peaks:
                self.reference_peaks = self.peak_results
            else:
                self.reference_peaks.update(self.peak_results)
        else:
            # Non-reference samples handling
            self.auto_select_peaks()
            self.fig.canvas.mpl_connect("key_press_event", self.on_key)
            self.fig.canvas.mpl_connect("button_press_event", self.on_click)
            plt.show(block=True)  # Blocks script until plot window is closed
        return self.peak_results, self.fig, self.reference_peaks, self.t_pressed

def run_peak_integrator(data, key, gi, pk_sns, smoothing_params, max_peaks_for_neighborhood, fp, gaussian_fit_mode, minimum_peak_amplitude=None, peak_prominence=0.001):
    # Setup data
    xdata = pd.Series(data['Samples'][key]['Raw Data'][data['Integration Metadata']['time_column']])
    ydata = pd.Series(data['Samples'][key]['Raw Data'][data['Integration Metadata']['signal_column']])

    # Subset to reference sample
    # --- Subset to global x-limits based on reference sample ---
    peak_times = list(data['Integration Metadata']['peak dictionary'].values())
    rt_buffer = 0.5  # 30 seconds = 0.5 minutes (you suggested 0.4, which is ~24s)

    xmin = min(peak_times) - rt_buffer
    xmax = max(peak_times) + rt_buffer
    mask = (xdata >= xmin) & (xdata <= xmax)

    xdata = xdata[mask].reset_index(drop=True)
    ydata = ydata[mask].reset_index(drop=True)

    ydata[ydata<0] = 0
    peak_timing = data['Integration Metadata']['peak dictionary'].values()
    data['Samples'][key]['Processed Data'] = {}

    # Match the HPLC/manual-FID background correction path: ASLS baseline,
    # clip negative corrected signal, then smooth for peak detection/fitting.
    base, min_peak_amp = hplc_style_baseline(xdata, ydata)
    y_bcorr = np.clip(ydata - base, 0, None)
    y_bcorr = smoother(pd.Series(y_bcorr, index=xdata.index), smoothing_params[0], smoothing_params[1])
    y_bcorr = pd.Series(y_bcorr, index=xdata.index)
    min_peak_amp = minimum_peak_amplitude if minimum_peak_amplitude is not None else min_peak_amp
    peak_indices, peak_properties = find_peaks(y_bcorr, height=min_peak_amp, prominence=peak_prominence)
    used_peaks = set()
    matched_indices = []
    presence_flags = []

    for pt in peak_timing:
        # Find candidate matches within tolerance
        distances = np.abs(xdata.iloc[peak_indices] - pt)
        candidates = [(idx, dist) for idx, dist in zip(peak_indices, distances) if dist <= 5/60]

        # Sort by closeness
        candidates.sort(key=lambda x: x[1])

        # Find the closest unused one
        selected = None
        for idx, dist in candidates:
            if idx not in used_peaks:
                selected = idx
                used_peaks.add(idx)
                break

        if selected is not None:
            matched_indices.append(selected)
            presence_flags.append(True)
        else:
            matched_indices.append(None)
            presence_flags.append(False)

    matched_indices = list(matched_indices)

    fig = plt.figure()
    plt.plot(xdata, y_bcorr, c= 'k', linewidth=1, linestyle='-', zorder=2)
    plt.title(str(key))
    valleys = find_valleys(y_bcorr, peak_indices)
    peak_labels = list(data['Integration Metadata']['peak dictionary'])
    for label, peak_idx in zip(peak_labels, matched_indices):
        if peak_idx is None:          # in case some peaks weren’t matched
            data['Samples'][key]['Processed Data'][label] = [np.nan]
            continue
        try:
            if gaussian_fit_mode in {"multi", "both", "asymmetric_or_multi"}:
                A, B, peak_neighborhood = find_peak_neighborhood_boundaries(
                    x=xdata, y_smooth=y_bcorr, peaks=peak_indices, valleys=valleys,
                    peak_idx=peak_idx, max_peaks=max_peaks_for_neighborhood,
                    peak_properties=peak_properties, gi=gi,
                    smoothing_params=smoothing_params, pk_sns=pk_sns)
                if not peak_neighborhood:
                    peak_neighborhood = [peak_idx]
            else:
                peak_neighborhood = [peak_idx]
            x_fit, y_fit_smooth, area_smooth, model_parameters = fit_gaussians(
                xdata, y_bcorr, peak_idx, peak_neighborhood,
                smoothing_params, pk_sns, gi=gi, mode=gaussian_fit_mode,
                valleys=valleys)
            plt.fill_between(x_fit, 0, y_fit_smooth, color=fit_color(model_parameters), alpha=0.5, zorder=1)
            x_peak_label = x_fit[np.argmax(y_fit_smooth)]
            add_peak_label(plt.gca(), label, x_fit, y_fit_smooth, xdata, y_bcorr)
            data['Samples'][key]['Processed Data'][label] = build_processed_peak_result(
                area_smooth, model_parameters, x_peak_label, x_fit, y_fit_smooth)
        except Exception as e:
            tqdm.write(f"[Warning] Failed to fit {label} in {key}: {e}")
            data['Samples'][key]['Processed Data'][label] = [np.nan]


    peak_times = list(data['Integration Metadata']['peak dictionary'].values())
    mean_val = np.mean(peak_times)
    xmin = min(peak_times) - mean_val * 0.1
    xmax = max(peak_times) + mean_val * 0.1

    # new y max
    mask = (xdata >= xmin) & (xdata <= xmax)
    y_max = ydata[mask].max()
    plt.xlim(xmin, xmax)
    plt.ylim(0, y_max+y_max*0.1)
    plt.ylabel(data['Integration Metadata']['signal_column'])
    plt.xlabel(data['Integration Metadata']['time_column'])
    plt.savefig(str(fp)+f"/{key}.png", dpi=300)
    plt.close()
    return data
