"""Onset determination: Sigma threshold and Poisson-CUSUM"""
import numpy as np
import pandas as pd
import pytest

from spearhead.vda import VDA, VDA_parameters
from spearhead.vda import analysis


def intensity_series(rise_at=30, spike_at=None, n_points=60):
    """5 min series with a noisy background and a rise to 10 at the rise_at point"""
    index = pd.date_range("2021-10-28 14:00", periods=n_points, freq="5min")
    values = np.where(np.arange(n_points) % 2 == 0, 1.0, 1.2)
    if rise_at is not None:
        values[rise_at:] = 10.0
    if spike_at is not None:
        values[spike_at] = 10.0
    return pd.Series(values, index=index)


def test_onset_detection_sigma_finds_rise():
    series = intensity_series(rise_at=30)
    onset, bg_start, bg_end, method_specific = VDA(VDA_parameters())._onset_detection_sigma(series, 3, 3, 0, 12)
    # the last point before the first point above the threshold
    assert onset == series.index[29]
    assert (bg_start, bg_end) == (series.index[0], series.index[12])
    assert method_specific["bg_level"] == pytest.approx(series.iloc[:13].mean())
    assert method_specific["threshold"] == pytest.approx(series.iloc[:13].mean() + 3 * series.iloc[:13].std())


def test_onset_detection_sigma_needs_consecutive_points():
    series = intensity_series(rise_at=40, spike_at=25)
    onset, *_ = VDA(VDA_parameters())._onset_detection_sigma(series, 3, 3, 0, 12)
    assert onset == series.index[39]


def test_onset_detection_sigma_without_onset():
    series = intensity_series(rise_at=None)
    onset, *_ = VDA(VDA_parameters())._onset_detection_sigma(series, 3, 3, 0, 12)
    assert onset is None


def test_onset_detection_sigma_rise_at_the_first_point():
    # without a point before the rise, the onset is its first point
    series = intensity_series(rise_at=None)
    series.iloc[:3] = 10.0
    onset, *_ = VDA(VDA_parameters())._onset_detection_sigma(series, 3, 3, 20, 40)
    assert onset == series.index[0]


def test_onset_detection_sigma_background_as_times():
    series = intensity_series(rise_at=30)
    v = VDA(VDA_parameters())
    by_index = v._onset_detection_sigma(series, 3, 3, 0, 12)
    by_time = v._onset_detection_sigma(series, 3, 3, series.index[0].to_pydatetime(), series.index[12].to_pydatetime())
    assert by_index[0] == by_time[0]
    assert by_index[3] == by_time[3]


def test_onset_detection_poisson_cusum_finds_rise():
    series = intensity_series(rise_at=30)
    onset, bg_start, bg_end, method_specific = VDA(VDA_parameters())._onset_detection_poisson_cusum(series, 30, 2, 0, 12)
    # pyonset (as SEPpy) gives the last point before the first alarm of the CUSUM function
    assert onset == series.index[29]
    assert (bg_start, bg_end) == (series.index[0], series.index[12])
    background = series.iloc[:13]
    assert method_specific["bg_level"] == pytest.approx(background.mean())
    assert method_specific["mu_d"] == pytest.approx(background.mean() + 2 * np.std(background))
    assert {"k", "h"} <= set(method_specific)


@pytest.mark.parametrize("freq, cusum_minutes, window", [
    ("5min", 30, 6),
    ("1min", 30, 30),
    ("1min", 15, 15),
    ("10min", 5, 1),
])
def test_onset_detection_poisson_cusum_window_from_cadence(monkeypatch, freq, cusum_minutes, window):
    """The CUSUM window in data points follows the cadence of the series (also without resampling)"""
    calls = []

    def onset_determination(ma_sigma, series, cusum_window, avg_end, sigma_multiplier=2):
        calls.append(cusum_window)
        return [0, 0, 0, 1, None, None, pd.NaT]

    monkeypatch.setattr(analysis, "onset_determination", onset_determination)
    series = pd.Series(1.0, index=pd.date_range("2021-10-28 14:00", periods=100, freq=freq))
    VDA(VDA_parameters())._onset_detection_poisson_cusum(series, cusum_minutes, 2, 0, 12)
    assert calls == [window]


def test_onset_detection_poisson_cusum_without_onset():
    series = intensity_series(rise_at=None)
    onset, *_ = VDA(VDA_parameters())._onset_detection_poisson_cusum(series, 30, 2, 0, 12)
    assert onset is None


def test_onset_detection_unknown_method():
    with pytest.raises(ValueError, match="Unknown onset method 'cusum'"):
        VDA(VDA_parameters())._onset_detection(intensity_series(), "cusum")
