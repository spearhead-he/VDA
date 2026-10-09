"""Unit tests that run offline: synthetic data, no particle data downloads and no SPICE kernels.

To run them locally, go to the base directory of the repository and run:
pytest tests/test_units.py
"""
import dataclasses
import os
import time
from datetime import datetime, timedelta

import astropy.constants as const
import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from spearhead.vda import OnsetSelection, VDA, VDA_parameters
from spearhead.vda import analysis, observers


def make_vda(**parameters) -> VDA:
    p = VDA_parameters()
    for name, value in parameters.items():
        setattr(p, name, value)
    return VDA(p)


def write_csv(tmp_path, name, text) -> str:
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# ---------------------------------------------------------------- onset detection

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


# ---------------------------------------------------------------- channel grouping

def test_group_channels_weighted_by_bin_width():
    index = pd.date_range("2021-10-28 14:00", periods=3, freq="5min")
    df = pd.DataFrame({"H_Flux_0": [1.0, 2.0, 3.0], "H_Flux_1": [3.0, 4.0, 5.0], "H_Flux_2": [10.0, 10.0, 10.0]}, index=index)
    grouped = VDA(VDA_parameters())._group_channels_de(
        df, [["H_Flux_0", "H_Flux_1"], ["H_Flux_2"]], [[1.0, 3.0], [2.0]], ["group 1", "group 2"]
    )
    # I = sum(I_n * dE_n) / sum(dE_n)
    assert grouped["group 1"].tolist() == pytest.approx([(1 * 1 + 3 * 3) / 4, (2 * 1 + 4 * 3) / 4, (3 * 1 + 5 * 3) / 4])
    assert grouped["group 2"].tolist() == pytest.approx([10.0, 10.0, 10.0])


# ---------------------------------------------------------------- events input

def test_date_ranges_without_file():
    v = make_vda(date_ranges=[(datetime(2021, 1, 1, 10), datetime(2021, 1, 1, 16)),
                              (datetime(2021, 2, 1, 8), datetime(2021, 2, 1, 12))])
    v.construct_times_df()
    assert v.df_times.index.tolist() == [1, 2]
    assert v.df_times.index.name == "Event No"
    assert v.df_times.loc[2, "Start Time"] == pd.Timestamp("2021-02-01 08:00")
    # default background window: 0 to 60 minutes after the start time
    assert v.bg_window(2) == (datetime(2021, 2, 1, 8), datetime(2021, 2, 1, 9))
    assert v.bg_window_source(2) == "default"


def test_date_ranges_start_after_end():
    v = make_vda(date_ranges=[(datetime(2021, 1, 2), datetime(2021, 1, 1))])
    with pytest.raises(ValueError, match="Event 1"):
        v.construct_times_df()


@pytest.mark.parametrize("header, row, expected_bg, source", [
    ("Event No,Start Time,End Time", "1,2021-05-22 19:45,2021-05-23 02:45",
     ("2021-05-22 19:45", "2021-05-22 20:45"), "default"),
    ("Event No,Start Time,BG End,End Time", "1,2021-05-22 19:45,2021-05-22 20:15,2021-05-23 02:45",
     ("2021-05-22 19:45", "2021-05-22 20:15"), "input file"),
    ("Event No, Start Time, BG Start, BG End, End Time",
     "1, 2021-05-22 19:45, 2021-05-22 20:00, 2021-05-22 20:30, 2021-05-23 02:45",
     ("2021-05-22 20:00", "2021-05-22 20:30"), "input file"),
    ("Event No,Start Time,BG Start,BG End,End Time", "1,2021-05-22 19:45,,,2021-05-23 02:45",
     ("2021-05-22 19:45", "2021-05-22 20:45"), "default"),
])
def test_datetime_ranges_file_layouts(tmp_path, header, row, expected_bg, source):
    path = write_csv(tmp_path, "events.csv", f"{header}\n{row}\n")
    v = make_vda(input_filepath=path)
    assert v.times_file_type(path) == "datetime ranges"
    v.construct_times_df()
    assert v.df_times.loc[1, "Start Time"] == pd.Timestamp("2021-05-22 19:45")
    assert v.df_times.loc[1, "End Time"] == pd.Timestamp("2021-05-23 02:45")
    assert v.bg_window(1) == tuple(pd.Timestamp(t).to_pydatetime() for t in expected_bg)
    assert v.bg_window_source(1) == source


def test_reference_times_file(tmp_path):
    path = write_csv(tmp_path, "events.csv", "Event,Reference Time\n7,2021-05-22 21:45\n")
    v = make_vda(input_filepath=path, bg_hours_prior=2, bg_hours_after=5)
    assert v.times_file_type(path) == "reference times"
    v.construct_times_df()
    event = v.df_times.loc[7]
    assert event["Start Time"] == pd.Timestamp("2021-05-22 19:45")
    assert event["End Time"] == pd.Timestamp("2021-05-23 02:45")
    assert event["Reference Time"] == pd.Timestamp("2021-05-22 21:45")


def test_events_file_with_unknown_columns(tmp_path):
    path = write_csv(tmp_path, "events.csv",
                     "Event No,BG Start,BG End,End Time\n"
                     "1,2021-05-22 19:45,2021-05-22 20:45,2021-05-23 02:45\n")
    v = make_vda(input_filepath=path)
    with pytest.raises(ValueError, match="Accepted columns"):
        v.construct_times_df()


# ---------------------------------------------------------------- background window

def test_background_window_changes(tmp_path):
    path = write_csv(tmp_path, "events.csv",
                     "Event No,Start Time,BG Start,BG End,End Time\n"
                     "1,2021-05-22 19:45,2021-05-22 20:00,2021-05-22 20:30,2021-05-23 02:45\n"
                     "2,2021-11-09 15:25,,,2021-11-09 22:25\n")
    v = make_vda(input_filepath=path)
    v.construct_times_df()

    v.set_default_bg_window(30, 90)
    # only the event without a window from the file changes
    assert v.bg_window(1) == (datetime(2021, 5, 22, 20, 0), datetime(2021, 5, 22, 20, 30))
    assert v.bg_window(2) == (datetime(2021, 11, 9, 15, 55), datetime(2021, 11, 9, 16, 55))

    v.set_bg_window(1, "2021-05-22 21:00", "2021-05-22 22:00")
    assert v.bg_window(1) == (datetime(2021, 5, 22, 21, 0), datetime(2021, 5, 22, 22, 0))
    assert v.bg_window_source(1) == "set"

    v.reset_bg_window(1)
    assert v.bg_window(1) == (datetime(2021, 5, 22, 20, 15), datetime(2021, 5, 22, 21, 15))
    assert v.bg_window_source(1) == "default"

    with pytest.raises(ValueError):
        v.set_bg_window(1, "2021-05-22 22:00", "2021-05-22 21:00")


def test_background_window_warnings():
    v = make_vda(date_ranges=[(datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 20))])
    v.construct_times_df()
    times = pd.date_range("2021-10-28 14:00", "2021-10-28 20:00", freq="5min")
    assert v.bg_window_warnings(1, times) == []
    v.set_default_bg_window(0, 5)
    assert any("2 data points" in w for w in v.bg_window_warnings(1, times))
    v.set_bg_window(1, "2021-10-28 19:00", "2021-10-28 21:00")
    assert any("outside its data range" in w for w in v.bg_window_warnings(1, times))


def test_save_times_round_trip(tmp_path):
    v = make_vda(date_ranges=[(datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 20)),
                              (datetime(2021, 11, 9, 15), datetime(2021, 11, 9, 22))])
    v.construct_times_df()
    v.set_bg_window(2, "2021-11-09 15:30", "2021-11-09 16:45")
    path = str(tmp_path / "saved.csv")
    v.save_times(path)

    saved = make_vda(input_filepath=path)
    saved.construct_times_df()
    columns = ["Start Time", "BG Start", "BG End", "End Time"]
    pd.testing.assert_frame_equal(saved.df_times[columns], v.df_times[columns], check_dtype=False)
    assert saved.bg_window_source(2) == "input file"


# ---------------------------------------------------------------- saved data

@pytest.mark.parametrize("saved_attrs, expected_frequency, message", [
    ({"resample_frequency": "10min"}, "10min", "Resample frequency set to '10min'"),
    ({"resample_frequency": "5min"}, "5min", None),
    ({}, "5min", "do not include their resample frequency"),
])
def test_loaded_data_resample_frequency(tmp_path, capsys, saved_attrs, expected_frequency, message):
    df = pd.DataFrame({"x": [1.0, 2.0]})
    df.attrs = dict(saved_attrs)
    df.to_pickle(tmp_path / "data.pkl")
    v = make_vda(resample_frequency="5min", load_data_filepath=str(tmp_path / "data.pkl"))
    v.construct_particles_df()
    assert v.parameters.resample_frequency == expected_frequency
    out = capsys.readouterr().out
    assert (message in out) if message else ("esample" not in out)


@pytest.mark.parametrize("saved_attrs", [{"observer": "solo"}, {}])
def test_loaded_data_observer(tmp_path, saved_attrs):
    # data saved without their observer are Solar Orbiter data
    df = pd.DataFrame({"x": [1.0, 2.0]})
    df.attrs = {"resample_frequency": "5min", **saved_attrs}
    df.to_pickle(tmp_path / "data.pkl")
    v = make_vda(load_data_filepath=str(tmp_path / "data.pkl"))
    v.construct_particles_df()


def test_loaded_data_of_another_observer(tmp_path):
    df = pd.DataFrame({"x": [1.0, 2.0]})
    df.attrs = {"observer": "other", "resample_frequency": "5min"}
    df.to_pickle(tmp_path / "data.pkl")
    v = make_vda(load_data_filepath=str(tmp_path / "data.pkl"))
    with pytest.raises(ValueError, match="data of the 'other' observer"):
        v.construct_particles_df()


# ---------------------------------------------------------------- observers

def epd_like_data(index, sensor="het"):
    """Data as returned by epd_load: protons with 3 channels and electrons with 2"""
    proton_flux, proton_bins = {"het": ("H_Flux", "H_Bins"), "ept": ("Ion_Flux", "Ion_Bins")}[sensor]

    def particle(flux_column, n_channels):
        columns = pd.MultiIndex.from_tuples(
            [(flux_column, f"{flux_column}_{c}") for c in range(n_channels)]
            + [("Uncertainty", f"Uncertainty_{c}") for c in range(n_channels)]
        )
        return pd.DataFrame(np.arange(len(index) * 2 * n_channels, dtype=float).reshape(len(index), -1),
                            index=index, columns=columns)

    energies = {
        f"{proton_bins}_Low_Energy": np.array([1.0, 2.0, 4.0]),
        f"{proton_bins}_Width": np.array([1.0, 2.0, 4.0]),
        "Electron_Bins_Low_Energy": np.array([0.5, 1.0]),
        "Electron_Bins_Width": np.array([0.5, 1.0]),
    }
    return particle(proton_flux, 3), particle("Electron_Flux", 2), energies


def test_solar_orbiter_load(monkeypatch):
    index = pd.date_range("2021-10-28 14:00", periods=4, freq="1min")
    monkeypatch.setattr(observers, "epd_load", lambda **kwargs: epd_like_data(index, kwargs["sensor"]))
    data = observers.SolarOrbiter().load("ept", index[0], index[-1], "sun", "path")
    assert list(data) == ["protons", "electrons"]
    assert data["protons"].flux.columns.tolist() == [0, 1, 2]
    assert data["protons"].flux[2].tolist() == epd_like_data(index, "ept")[0][("Ion_Flux", "Ion_Flux_2")].tolist()
    assert data["electrons"].energies.loc[1].tolist() == [1.0, 1.0]


def test_data_from_the_observer(monkeypatch):
    # the observer data become the columns of the analysis, e.g. H_Flux_2 for proton channel 2
    index = pd.date_range("2021-10-28 14:00", periods=20, freq="1min")
    monkeypatch.setattr(observers, "epd_load", lambda **kwargs: epd_like_data(index, kwargs["sensor"]))
    v = make_vda(date_ranges=[(index[0], index[-1])], resample_frequency="5min", viewings=["sun", "north"],
                 channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}}})
    v.construct_times_df()
    v.construct_energies_df()
    assert v.df_energies.loc[("het", "H_Flux_2")].tolist() == [4.0, 4.0, 8.0]
    assert v.df_energies.loc[("ept", "Electron_Flux_1"), "High Energy"] == 2.0
    v.construct_particles_df()
    assert v.df_data.attrs["observer"] == "solo"
    assert v.df_data.columns.tolist() == [("het", "protons", viewing, "H_Flux", f"H_Flux_{c}")
                                          for viewing in ("sun", "north") for c in range(3)]
    assert len(v.df_data) == 4


def stereo_like_data(instrument, startdate, enddate, spacecraft, sept_species=None, sept_viewing=None, path=None):
    """Data as returned by stereo_load: HET with 2 proton and 2 electron channels, SEPT with channels 2 to 4"""
    index = pd.date_range("2021-10-28 14:00", periods=20, freq="1min")

    def channels_dict_df(channels):
        lower = np.array([1.0, 2.0, 4.0])[:len(channels)]
        return pd.DataFrame({"lower_E": lower, "upper_E": 2 * lower, "DE": lower}, index=channels)

    if instrument == "HET":
        columns = ["Electron_Flux_0", "Electron_Flux_1", "Electron_Sigma_0", "Electron_Sigma_1",
                   "Proton_Flux_0", "Proton_Flux_1", "Proton_Sigma_0", "Proton_Sigma_1"]
        meta = {"channels_dict_df_e": channels_dict_df([0, 1]), "channels_dict_df_p": channels_dict_df([0, 1])}
    else:
        assert sept_viewing in ("sun", "asun", "north", "south")
        columns = ["ch_2", "ch_3", "ch_4", "err_ch_2", "err_ch_3", "err_ch_4"]
        meta = {f"channels_dict_df_{sept_species}": channels_dict_df([2, 3, 4])}
    flux = pd.DataFrame(np.ones((len(index), len(columns))), index=index, columns=columns)
    return flux, meta


def test_stereo_a_load(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "stereo_load", stereo_like_data)
    sta = observers.StereoA()
    het = sta.load("het", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path))
    assert list(het) == ["protons", "electrons"]
    assert het["protons"].flux.columns.tolist() == [0, 1]
    assert het["electrons"].energies.loc[1].tolist() == [2.0, 2.0]
    sept = sta.load("sept", "2021-10-28 14:00", "2021-10-28 15:00", "north", str(tmp_path), particles=["electrons"])
    assert list(sept) == ["electrons"]
    assert sept["electrons"].flux.columns.tolist() == [2, 3, 4]
    assert sept["electrons"].energies.loc[4].tolist() == [4.0, 4.0]


def test_stereo_a_without_data(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "stereo_load", lambda *args, **kwargs: ([], []))
    with pytest.raises(ValueError, match="No STEREO-A HET data"):
        observers.StereoA().load("het", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path))


def test_viewings_of_each_sensor(monkeypatch):
    # STEREO-A HET has only the omni viewing, SEPT has no omni viewing
    monkeypatch.setattr(observers, "stereo_load", stereo_like_data)
    v = make_vda(observer="sta", viewings=["north", "omni", "sun"],
                 channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [0, 1]}},
                                 "electrons": {"SEPT/electrons Channel 1": {"sensor": "sept", "channels": [2, 3]}}})
    assert v.sensor_viewings("het") == ["omni"]
    assert v.sensor_viewings("sept") == ["north", "sun"]
    v.construct_times_df()
    v.construct_energies_df()
    v.construct_particles_df()
    assert sorted({c[:3] for c in v.df_data.columns}) == [
        ("het", "protons", "omni"), ("sept", "electrons", "north"), ("sept", "electrons", "sun")]
    v.group_energy_channels()
    v.construct_energy_channels_characteristics()
    assert v.df_channels_chars.index.get_level_values("sensor").tolist() == ["het", "sept"]


def test_sensor_without_selected_viewings():
    v = make_vda(observer="sta", viewings=["sun"])
    with pytest.raises(ValueError, match="None of the selected viewings \\(sun\\) is a viewing of STEREO-A HET"):
        v.check_viewings()


def test_observer_change_resets_viewings_and_channel_groups():
    p = VDA_parameters(viewings=["north"])
    p.observer = "solo"
    assert p.viewings == ["north"]
    p.observer = "sta"
    assert p.viewings == ["sun", "omni"]
    assert p.sensors_particles == {"het": ["protons", "electrons"]}
    assert [g["channels"] for g in p.channel_groups["protons"].values()] == [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10]]
    assert [g["channels"] for g in p.channel_groups["electrons"].values()] == [[0, 1, 2]]
    assert p.AVAILABLE_CHANNELS["sept"]["electrons"] == tuple(range(2, 17))


@pytest.mark.parametrize("n_channels, particle, expected", [
    (6, "protons", [[0, 1, 2], [3, 4, 5]]),
    (8, "protons", [[0, 1, 2], [3, 4, 5], [6, 7]]),
    (7, "protons", [[0, 1, 2], [3, 4, 5, 6]]),
    (4, "electrons", [[0, 1], [2, 3]]),
    (5, "electrons", [[0, 1], [2, 3, 4]]),
    (1, "electrons", [[0]]),
])
def test_default_groups(n_channels, particle, expected):
    # a single channel left at the end joins the previous group
    assert observers.default_groups(tuple(range(n_channels)), particle) == expected


def test_unknown_observer():
    with pytest.raises(ValueError, match="Unknown observer 'other'"):
        VDA_parameters(observer="other")
    p = VDA_parameters()
    with pytest.raises(ValueError, match="Unknown observer 'other'"):
        p.observer = "other"


# ---------------------------------------------------------------- onset selection

def vda_with_onset_options(viewings):
    """Two channels: onsets in sun and north for the first, only in north for the second"""
    v = make_vda(viewings=viewings)
    rows = [(1, "het", "protons", "H_Flux", "HET/protons Channel 1", "sun"),
            (1, "het", "protons", "H_Flux", "HET/protons Channel 1", "north"),
            (1, "het", "protons", "H_Flux", "HET/protons Channel 2", "north")]
    v.df_options = pd.DataFrame(
        {"Onset Time": [pd.Timestamp("2021-10-28 16:00")] * len(rows)},
        index=pd.MultiIndex.from_tuples(rows, names=["Event No", "sensor", "particle", "prefix", "channels", "viewing"]),
    )
    return v


@pytest.mark.parametrize("viewings, expected", [
    (["sun", "north"], ["sun", "north"]),
    (["north", "sun"], ["north", "north"]),
])
def test_use_all_follows_viewings_order(viewings, expected):
    v = vda_with_onset_options(viewings)
    v.select_onsets()
    assert v.parameters.selected_onsets["Viewing"].tolist() == expected


# ---------------------------------------------------------------- VDA fit

class FixedDistance:
    distance = 0.8 * u.AU


def vda_with_onsets(onsets):
    """One event with onsets on 4 channels of inverse beta 1 to 4"""
    channels = [f"HET/protons Channel {i}" for i in range(1, len(onsets) + 1)]
    v = make_vda(date_ranges=[(datetime(2021, 10, 31, 2), datetime(2021, 10, 31, 6))])
    v.construct_times_df()
    v.df_channels_chars = pd.DataFrame(
        {"Inverse Beta": [float(i) for i in range(1, len(onsets) + 1)]},
        index=pd.MultiIndex.from_tuples([("het", "protons", c) for c in channels], names=["sensor", "particle", "channel"]),
    )
    v.df_onsets_existing = pd.DataFrame(
        {"Onset Time": [pd.Timestamp(t) for t in onsets]},
        index=pd.MultiIndex.from_tuples([(1, "het", "protons", "sun", "H_Flux", c) for c in channels]),
    )
    v.df_options = v.df_onsets_existing.reorder_levels([0, 1, 2, 4, 5, 3])
    v.parameters.selected_onsets = pd.DataFrame(
        {"Viewing": ["sun"] * len(channels)},
        index=pd.MultiIndex.from_tuples([(1, "het", "protons", "H_Flux", c) for c in channels]),
    )
    return v


@pytest.fixture
def fixed_spacecraft_distance(monkeypatch):
    monkeypatch.setattr(observers.SolarOrbiter, "distance", lambda self, time: FixedDistance.distance)
    return (FixedDistance.distance / const.c).to(u.s).value


def check_exact_line_fit(extra_seconds):
    # onsets exactly on a line: 1200 s per inverse beta, crossing inverse beta 0 at 03:00
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00", "2021-10-31 04:20"])
    v.compute_vda()
    result = v.results.loc[1]
    expected_release = datetime(2021, 10, 31, 3, 0) + timedelta(seconds=extra_seconds)
    assert result["Release Time"] == expected_release.strftime("%Y-%m-%d %H:%M:%S")
    assert pd.Timedelta(result["Extra Time"]).total_seconds() == pytest.approx(extra_seconds, abs=1e-5)
    assert result["APL"] == pytest.approx(1200 / extra_seconds)
    assert pd.Timedelta(result["Release Time Error"]).total_seconds() == pytest.approx(0, abs=1e-3)


def test_vda_fit_exact_line(fixed_spacecraft_distance):
    check_exact_line_fit(fixed_spacecraft_distance)


@pytest.mark.skipif(not hasattr(time, "tzset"), reason="time.tzset is not available on this platform")
def test_vda_fit_independent_of_local_timezone(fixed_spacecraft_distance):
    # the onsets are in the repeated hour of the 31 October 2021 daylight saving change in Athens
    tz = os.environ.get("TZ")
    try:
        os.environ["TZ"] = "Europe/Athens"
        time.tzset()
        check_exact_line_fit(fixed_spacecraft_distance)
    finally:
        if tz is None:
            os.environ.pop("TZ", None)
        else:
            os.environ["TZ"] = tz
        time.tzset()


def test_vda_fit_needs_two_points(fixed_spacecraft_distance, capsys):
    v = vda_with_onsets(["2021-10-31 03:20"])
    v.compute_vda()
    assert "Not enough onset points in event 1" in capsys.readouterr().out
    assert pd.isna(v.results.loc[1, "APL"])


# ---------------------------------------------------------------- parameters

@pytest.mark.parametrize("name, version", [
    ("input_type", "v0.3.0"), ("date_range_filepath", "v0.3.0"), ("reference_times_filepath", "v0.3.0"),
    ("date_start", "v0.3.0"), ("date_end", "v0.3.0"),
    ("viewings_tt", "v0.4.0"), ("load_data", "v0.4.0"), ("save_data", "v0.4.0"),
    ("default_channel_groups", "v0.5.0"),
])
def test_removed_parameters(name, version):
    with pytest.raises(AttributeError, match=f"removed in {version}"):
        setattr(VDA_parameters(), name, None)


def test_unknown_parameter():
    with pytest.raises(AttributeError, match="no parameter 'input_filepth'"):
        VDA_parameters().input_filepth = "events.csv"


def test_onset_selection_values():
    p = VDA_parameters()
    p.onset_selection = 1
    assert p.onset_selection is OnsetSelection.INTERACTIVE
    with pytest.raises(ValueError, match="Unknown onset selection 2"):
        p.onset_selection = 2


def test_parameters_instances_are_independent():
    p = VDA_parameters()
    q = dataclasses.replace(p, input_filepath="events.csv")
    assert p.input_filepath == "" and q.input_filepath == "events.csv"
    # new parameter sets do not share their lists and dictionaries
    other = VDA_parameters()
    other.date_ranges.append((datetime(2021, 1, 1), datetime(2021, 1, 2)))
    other.channel_groups["protons"]["HET/protons Channel 1"]["channels"].append(4)
    assert len(p.date_ranges) == 1 and p.channel_groups["protons"]["HET/protons Channel 1"]["channels"] == [1, 2, 3]
    assert VDA_parameters() == VDA_parameters()


def test_onset_method_parameters_follow_the_method():
    p = VDA_parameters()
    assert p.onset_method_parameters == {"s": 3, "n": 3}
    p.onset_method_parameters["s"] = 4
    p.onset_method = "sigma"
    assert p.onset_method_parameters == {"s": 4, "n": 3}
    p.onset_method = "poisson_cusum"
    assert p.onset_method_parameters == {"cusum_minutes": 30, "sigma_multiplier": 2}
    assert VDA_parameters(onset_method="poisson_cusum").onset_method_parameters == {"cusum_minutes": 30, "sigma_multiplier": 2}
    with pytest.raises(ValueError, match="Unknown onset method 'cusum'"):
        p.onset_method = "cusum"


def test_onset_method_parameters_of_another_method():
    v = make_vda()
    v.parameters.onset_method_parameters = {"cusum_minutes": 30}
    with pytest.raises(ValueError, match="not parameters of the sigma onset method"):
        v.calculate_onsets()


def test_default_channel_groups():
    p = VDA_parameters()
    assert list(p.channel_groups) == ["protons", "electrons"]
    assert len(p.channel_groups["protons"]) == 8
    assert p.channel_groups["electrons"] == {
        "HET/electrons Channel 1": {"sensor": "het", "channels": [0, 1]},
        "HET/electrons Channel 2": {"sensor": "het", "channels": [2, 3]},
    }
    assert p.sensors_particles == {"het": ["protons", "electrons"]}
    assert p.viewings == ["sun"]
