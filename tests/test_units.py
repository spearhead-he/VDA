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

import vda as vda_module
from vda import VDA
from vda_tool_configuration import OnsetSelection, VDA_parameters
from vda_views import VDA_nb_displayer


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
    assert onset == series.index[30]
    assert (bg_start, bg_end) == (series.index[0], series.index[12])
    assert method_specific["bg_level"] == pytest.approx(series.iloc[:13].mean())
    assert method_specific["threshold"] == pytest.approx(series.iloc[:13].mean() + 3 * series.iloc[:13].std())


def test_onset_detection_sigma_needs_consecutive_points():
    series = intensity_series(rise_at=40, spike_at=25)
    onset, *_ = VDA(VDA_parameters())._onset_detection_sigma(series, 3, 3, 0, 12)
    assert onset == series.index[40]


def test_onset_detection_sigma_without_onset():
    series = intensity_series(rise_at=None)
    onset, *_ = VDA(VDA_parameters())._onset_detection_sigma(series, 3, 3, 0, 12)
    assert onset is None


def test_onset_detection_sigma_background_as_times():
    series = intensity_series(rise_at=30)
    v = VDA(VDA_parameters())
    by_index = v._onset_detection_sigma(series, 3, 3, 0, 12)
    by_time = v._onset_detection_sigma(series, 3, 3, series.index[0].to_pydatetime(), series.index[12].to_pydatetime())
    assert by_index[0] == by_time[0]
    assert by_index[3] == by_time[3]


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
    assert v._bg_window(2) == (datetime(2021, 2, 1, 8), datetime(2021, 2, 1, 9))
    assert v._bg_window_source(2) == "default"


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
    assert v._bg_window(1) == tuple(pd.Timestamp(t).to_pydatetime() for t in expected_bg)
    assert v._bg_window_source(1) == source


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
    assert v._bg_window(1) == (datetime(2021, 5, 22, 20, 0), datetime(2021, 5, 22, 20, 30))
    assert v._bg_window(2) == (datetime(2021, 11, 9, 15, 55), datetime(2021, 11, 9, 16, 55))

    v.set_bg_window(1, "2021-05-22 21:00", "2021-05-22 22:00")
    assert v._bg_window(1) == (datetime(2021, 5, 22, 21, 0), datetime(2021, 5, 22, 22, 0))
    assert v._bg_window_source(1) == "set"

    v.reset_bg_window(1)
    assert v._bg_window(1) == (datetime(2021, 5, 22, 20, 15), datetime(2021, 5, 22, 21, 15))
    assert v._bg_window_source(1) == "default"

    with pytest.raises(ValueError):
        v.set_bg_window(1, "2021-05-22 22:00", "2021-05-22 21:00")


def test_background_window_warnings():
    v = make_vda(date_ranges=[(datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 20))])
    v.construct_times_df()
    times = pd.date_range("2021-10-28 14:00", "2021-10-28 20:00", freq="5min")
    assert v._bg_window_warnings(1, times) == []
    v.set_default_bg_window(0, 5)
    assert any("2 data points" in w for w in v._bg_window_warnings(1, times))
    v.set_bg_window(1, "2021-10-28 19:00", "2021-10-28 21:00")
    assert any("outside its data range" in w for w in v._bg_window_warnings(1, times))


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
    assert saved._bg_window_source(2) == "input file"


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
    VDA_nb_displayer(v).select_onsets()
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
    monkeypatch.setattr(vda_module.spice, "get_body", lambda *args, **kwargs: FixedDistance())
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

@pytest.mark.parametrize("name", ["input_type", "date_range_filepath", "reference_times_filepath",
                                  "date_start", "date_end", "viewings_tt"])
def test_removed_parameters(name):
    with pytest.raises(AttributeError, match="removed in v0.3.0"):
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
    other.channel_groups["protons"] = {}
    assert len(p.date_ranges) == 1 and p.channel_groups == {}
    assert VDA_parameters() == VDA_parameters()
