"""Events input and background windows"""
from datetime import datetime

import pandas as pd
import pytest

from tests.helpers import make_vda, write_csv


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


def test_events_file_not_comma_separated(tmp_path):
    path = write_csv(tmp_path, "events.csv",
                     "Event No;Start Time;End Time\n"
                     "1;2021-05-22 19:45;2021-05-23 02:45\n")
    v = make_vda(input_filepath=path)
    with pytest.raises(ValueError, match="The columns of .*events.csv are separated by ';': separate them with commas"):
        v.construct_times_df()


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


def test_events_file_not_found(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    v = make_vda(input_filepath="events.csv")
    with pytest.raises(FileNotFoundError, match="Events file not found: events.csv .* relative to the current folder"):
        v.construct_times_df()


def test_events_file_with_invalid_datetimes(tmp_path):
    path = write_csv(tmp_path, "events.csv", "Event No,Start Time,End Time\n"
                                             "1,2021-10-28 14:00,2021-10-28 20:00\n"
                                             "2,2021-11-09 15:00,yesterday\n")
    v = make_vda(input_filepath=path)
    with pytest.raises(ValueError, match=r"Invalid datetimes in the column End Time of .*events.csv, in the rows "
                                         r"2 \('yesterday'\)"):
        v.construct_times_df()


def test_events_file_with_empty_background_cells(tmp_path):
    # empty cells are not invalid datetimes: they use the default background window
    path = write_csv(tmp_path, "events.csv", "Event No,Start Time,BG Start,BG End,End Time\n"
                                             "1,2021-10-28 14:00,,,2021-10-28 20:00\n")
    v = make_vda(input_filepath=path)
    v.construct_times_df()
    assert v.bg_window_source(1) == "default"


def test_events_file_start_after_end(tmp_path):
    path = write_csv(tmp_path, "events.csv", "Event No,Start Time,End Time\n1,2021-10-28 20:00,2021-10-28 14:00\n")
    v = make_vda(input_filepath=path)
    with pytest.raises(ValueError, match="Event 1 of .*events.csv: the datetime range start .* must be before its end"):
        v.construct_times_df()
