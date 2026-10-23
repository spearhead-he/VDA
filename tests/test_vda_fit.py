"""VDA fit and plot"""
import os
import time
from datetime import datetime, timedelta

import astropy.constants as const
import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import PathCollection

from spearhead.vda import analysis, observers, views
from tests.helpers import legend_texts, make_vda, reference_legend, rising_epd_data, table_rows


class FixedDistance:
    distance = 0.8 * u.AU


def vda_with_onsets(onsets, particles=None):
    """One event with onsets on channels of inverse beta 1, 2, ..., of HET protons or of the given particles"""
    particles = particles or ["protons"] * len(onsets)
    channels = [f"HET/{p} Channel {i}" for i, p in enumerate(particles, start=1)]
    prefixes = ["H_Flux" if p == "protons" else "Electron_Flux" for p in particles]
    v = make_vda(date_ranges=[(datetime(2021, 10, 31, 2), datetime(2021, 10, 31, 6))])
    v.construct_times_df()
    v.df_channels_chars = pd.DataFrame(
        {"Inverse Beta": [float(i) for i in range(1, len(onsets) + 1)]},
        index=pd.MultiIndex.from_tuples(
            [("het", p, c) for p, c in zip(particles, channels)], names=["sensor", "particle", "channel"]
        ),
    )
    v.df_onsets_existing = pd.DataFrame(
        {"Onset Time": [pd.Timestamp(t) for t in onsets]},
        index=pd.MultiIndex.from_tuples([(1, "het", p, "sun", f, c) for p, f, c in zip(particles, prefixes, channels)]),
    )
    v.df_options = v.df_onsets_existing.reorder_levels([0, 1, 2, 4, 5, 3])
    v.parameters.selected_onsets = pd.DataFrame(
        {"Viewing": ["sun"] * len(channels)},
        index=pd.MultiIndex.from_tuples([(1, "het", p, f, c) for p, f, c in zip(particles, prefixes, channels)]),
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
    # the path length of the slope, in AU, whatever the distance of the observer
    assert result["APL"] == pytest.approx(1200 / (u.AU / const.c).to(u.s).value)
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
    assert pd.isna(v.results.loc[1, "APL"]) and v.vda_fits == {}
    v.print_results(1)
    assert "No results (not enough onset points: they are needed at two different energies at least)" \
        in capsys.readouterr().out


def test_vda_points_of_protons_and_electrons(fixed_spacecraft_distance):
    v = vda_with_onsets(
        ["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00", "2021-10-31 04:20"],
        particles=["protons", "electrons", "protons", "electrons"],
    )
    v.compute_vda()
    assert list(v.vda_fits[1]["particles"]) == ["protons", "electrons", "protons", "electrons"]
    # the date of the title
    times = pd.date_range("2021-10-31 02:00", periods=2, freq="min")
    v.df_grouped = pd.DataFrame(index=pd.MultiIndex.from_product([[1], times]))
    fig = views.plot_vda(v, 1)
    # the title names the spacecraft
    assert fig.get_suptitle() == "Event 1 (2021-10-31) | Solar Orbiter"
    # the points in black, the release time at the Sun as an open circle at inverse beta 0
    points = [c for c in fig.axes[0].collections if len(c.get_offsets()) and c.get_offsets()[0][0] != 0]
    assert [len(c.get_offsets()) for c in points] == [2, 2]
    assert [c.get_paths()[0] for c in points] != [points[0].get_paths()[0]] * 2
    assert legend_texts(fig.axes[0])[:2] == ["Protons", "Electrons"]
    plt.close(fig)


def test_vda_plot_reference_times(fixed_spacecraft_distance):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00"])
    v.compute_vda()
    times = pd.date_range("2021-10-31 02:00", periods=2, freq="min")
    v.df_grouped = pd.DataFrame(index=pd.MultiIndex.from_product([[1], times]))
    references = [("SXR peak", datetime(2021, 10, 31, 2, 50)), ("type II", pd.Timestamp("2021-10-31 03:00"), 0.5)]
    ax = views.plot_vda(v, 1, references=references).axes[0]
    # the reference times at inverse beta 0 at the Sun, with different markers, and their times in the legend
    labels = ["SXR peak: 02:41:44 (02:50:00 at Earth)", "type II: 02:55:50 (03:00:00 at 0.50 AU)"]
    points = {c.get_label(): c for c in ax.collections if c.get_label() in labels}
    assert list(points) == labels
    assert all(c.get_offsets()[0][0] == 0 for c in points.values())
    type_ii = mdates.num2date(points[labels[1]].get_offsets()[0][1]).replace(tzinfo=None)
    assert abs(type_ii - datetime(2021, 10, 31, 2, 55, 50, 497611)) < timedelta(milliseconds=1)
    assert points[labels[0]].get_paths()[0] != points[labels[1]].get_paths()[0]
    assert ax.get_xlim()[0] == 0
    assert [t.get_text() for t in reference_legend(ax).get_texts()] == labels
    # without reference times, no legend of them
    assert reference_legend(views.plot_vda(v, 1).axes[0]) is None
    plt.close("all")


def test_vda_plot_from_the_release_at_the_sun(fixed_spacecraft_distance):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00"])
    v.compute_vda()
    v.df_grouped = pd.DataFrame(index=pd.MultiIndex.from_product([[1], pd.date_range("2021-10-31 02:00", periods=2)]))
    ax = views.plot_vda(v, 1).axes[0]
    # the fit line from inverse beta 0, where an open circle marks the release time at the Sun
    fit_line = ax.get_lines()[0]
    assert fit_line.get_xdata()[0] == 0 and ax.get_xlim()[0] == 0
    at_sun = analysis.release_times(v.vda_fits[1])[0]
    release = [c for c in ax.collections if isinstance(c, PathCollection) and c.get_offsets()[0][0] == 0][0]
    for plotted in (mdates.num2date(release.get_offsets()[0][1]).replace(tzinfo=None), fit_line.get_ydata()[0]):
        assert abs(plotted - at_sun.to_pydatetime()) < timedelta(seconds=1)
    plt.close("all")


def test_release_times_at_the_sun_and_at_1_au(fixed_spacecraft_distance):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00", "2021-10-31 04:21"])
    v.compute_vda()
    fit = v.vda_fits[1]
    at_sun, at_au = analysis.release_times(fit)
    # to the microsecond
    assert abs(at_sun - (pd.Timestamp(fit["t0"]) + pd.Timedelta(seconds=fit["b"]))) <= pd.Timedelta(microseconds=1)
    assert at_sun == at_sun.round("us")
    # the light travel time over 1 AU, 8 min 19 s
    assert (at_au - at_sun).total_seconds() == pytest.approx(499.005, abs=1e-3)
    # a fit through 03:00:00 at inverse beta 0 is shown as 03:00:00, not 02:59:59
    assert analysis.release_times({"t0": datetime(2021, 10, 31, 2), "b": 3599.9999999})[0] == pd.Timestamp("2021-10-31 03:00")
    # the printed release time is the time at the Sun plus the light travel time to the spacecraft
    res = v.results.loc[1]
    assert pd.Timestamp(res["Release Time"]) - res["Extra Time"] == pytest.approx(at_sun, abs=pd.Timedelta(seconds=1))


def test_vda_fit_of_two_points_has_no_errors(fixed_spacecraft_distance, capsys):
    # the line goes through both points: the errors cannot be estimated
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40"])
    v.compute_vda()
    result = v.results.loc[1]
    assert result["Release Time"] == (datetime(2021, 10, 31, 3) + timedelta(seconds=fixed_spacecraft_distance)
                                      ).strftime("%Y-%m-%d %H:%M:%S")
    assert pd.isna(result["Release Time Error"]) and pd.isna(result["APL Error"])
    v.print_results(1)
    out = capsys.readouterr().out
    assert "    Release at the Sun : 2021-10-31 03:00:00 ± n/a\n" in out and " ± n/a AU\n" in out
    v.df_grouped = pd.DataFrame(index=pd.MultiIndex.from_product([[1], pd.date_range("2021-10-31 02:00", periods=3)]))
    fig = views.plot_vda(v, 1)
    # no band of the error
    assert "2σ band of the fit" not in fig.axes[0].get_legend_handles_labels()[1]
    plt.close(fig)


def test_vda_fit_needs_two_energies(fixed_spacecraft_distance):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40"])
    v.df_channels_chars["Inverse Beta"] = 2.0
    v.compute_vda()
    assert pd.isna(v.results.loc[1, "APL"]) and v.vda_fits == {}


@pytest.mark.parametrize("onsets, warnings", [
    # later onsets at the higher energies
    (["2021-10-31 04:20", "2021-10-31 04:00", "2021-10-31 03:40"], ["the APL is not positive"]),
    # 1 minute per inverse beta: an APL of 0.12 AU, at 0.8 AU from the Sun
    (["2021-10-31 03:01", "2021-10-31 03:02", "2021-10-31 03:03", "2021-10-31 03:04"],
     ["the APL (0.12 AU) is shorter than the distance of the observer from the Sun (0.80 AU)"]),
    # the line crosses inverse beta 0 at 03:50, after the earliest onset (03:40), with an APL of 0.72 AU
    (["2021-10-31 04:00", "2021-10-31 04:10", "2021-10-31 03:40", "2021-10-31 04:30"],
     ["the APL (0.72 AU) is shorter than the distance", "the release at the Sun is after the earliest onset"]),
])
def test_vda_fit_warnings(fixed_spacecraft_distance, capsys, onsets, warnings):
    v = vda_with_onsets(onsets)
    v.compute_vda()
    assert len(v.vda_fits[1]["warnings"]) == len(warnings)
    assert all(w.startswith(expected) for w, expected in zip(v.vda_fits[1]["warnings"], warnings))
    v.print_results(1)
    out = capsys.readouterr().out
    assert all(f"    Warning: {expected}" in out for expected in warnings)


def test_vda_fit_apl_of_the_distance_without_warnings():
    # a path length of the distance of the observer from the Sun (light travel time 400 s)
    assert analysis.fit_warnings(400.0, 0.0, np.array([400.0, 800.0]), 400.0) == []


def test_vda_fit_without_warnings(fixed_spacecraft_distance):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00", "2021-10-31 04:20"])
    v.compute_vda()
    assert v.vda_fits[1]["warnings"] == []


def test_printed_release_time_with_the_light_travel_time(fixed_spacecraft_distance, capsys):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00", "2021-10-31 04:21"])
    v.compute_vda()
    v.print_results(1)
    out = capsys.readouterr().out
    res = v.results.loc[1]
    # the release time at the Sun and its error first, then the times its light reaches 1 AU and the spacecraft
    assert out.splitlines()[1:5] == [
        f"    Release at the Sun : {res['Release Time (Sun)']} ± {analysis.format_timedelta(res['Release Time Error'])}",
        f"    At 1 AU            : {res['Release Time (1 AU)']}",
        f"    At the spacecraft  : {res['Release Time']} (+0:06:39 light travel time from the Sun)",
        f"    APL                : {res['APL']:.2f} ± {res['APL Error']:.2f} AU",
    ]
    assert " ± 0:00:" in out.splitlines()[1]
    v.df_grouped = pd.DataFrame(index=pd.MultiIndex.from_product([[1], pd.date_range("2021-10-31 02:00", periods=3)]))
    fig = views.plot_vda(v, 1)
    assert "2σ band" in legend_texts(fig.axes[0])
    assert table_rows(fig)[1][3].startswith("± 0:00:")
    plt.close(fig)


def test_vda_results_of_the_events_of_the_last_run(fixed_spacecraft_distance):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00"])
    v.results.loc[2] = v.results.columns.map(lambda c: None)
    v.compute_vda()
    assert v.results.index.tolist() == [1]


def two_events_vda(monkeypatch):
    """Two events of the synthetic data, with one grouped channel of protons and one of electrons"""
    monkeypatch.setattr(observers, "epd_load", rising_epd_data)
    monkeypatch.setattr(analysis.VDA, "define_spacecraft_parameters", lambda self: None)
    monkeypatch.setattr(observers.SolarOrbiter, "distance", lambda self, time: FixedDistance.distance)
    return make_vda(
        date_ranges=[(datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 20)),
                     (datetime(2021, 11, 9, 15), datetime(2021, 11, 9, 21))],
        channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}},
                        "electrons": {"HET/electrons Channel 1": {"sensor": "het", "channels": [0, 1]}}},
    )


def test_run_gives_the_results_of_the_steps(monkeypatch):
    v_steps = two_events_vda(monkeypatch)
    for step in ("construct_times_df", "construct_energies_df", "construct_particles_df", "group_energy_channels",
                 "calculate_onsets", "clean_onsets", "construct_options_df", "select_onsets",
                 "construct_energy_channels_characteristics", "define_spacecraft_parameters", "compute_vda"):
        getattr(v_steps, step)()
    v = two_events_vda(monkeypatch)
    results = v.run()
    assert results is v.results
    pd.testing.assert_frame_equal(results, v_steps.results)
    # the tables of the steps stay available
    for table in ("df_times", "df_energies", "df_data", "df_grouped", "df_onsets_existing", "df_options",
                  "df_channels_chars"):
        pd.testing.assert_frame_equal(getattr(v, table), getattr(v_steps, table))
    assert v.vda_fits.keys() == v_steps.vda_fits.keys() == {1, 2}


def test_run_with_background_windows(monkeypatch):
    v = two_events_vda(monkeypatch)
    v.run(bg_windows={2: ("2021-11-09 15:30", "2021-11-09 16:30")})
    assert v.bg_window(2) == (datetime(2021, 11, 9, 15, 30), datetime(2021, 11, 9, 16, 30))
    assert (v.bg_window_source(1), v.bg_window_source(2)) == ("default", "set")
    background = v.df_onsets_existing.xs(2, level=0)["Background Start"]
    assert (background == pd.Timestamp("2021-11-09 15:30")).all()


def test_background_window_of_an_unknown_event(monkeypatch):
    v = two_events_vda(monkeypatch)
    with pytest.raises(ValueError, match="No event 3: the events are 1, 2"):
        v.run(bg_windows={3: ("2021-11-09 15:30", "2021-11-09 16:30")})


def test_vda_steps_in_order():
    v = make_vda()
    with pytest.raises(RuntimeError, match=r"Run the onsets step \(tool.run_onsets\(\) in the notebook\) first"):
        v.compute_vda()
    with pytest.raises(RuntimeError, match=r"Run the data step \(tool.run_data\(\) in the notebook\) first"):
        v.calculate_onsets()
    with pytest.raises(RuntimeError, match="Run the data step"):
        v.construct_energy_channels_characteristics()


def test_release_time_text():
    res = pd.Series({"Release Time": "2021-10-28 15:26:11", "Release Time Error": timedelta(minutes=3, seconds=38.9),
                     "Extra Time": timedelta(minutes=6, seconds=40.9), "APL": 1.76, "APL Error": 0.12})
    assert analysis.format_release_time(res) == "2021-10-28 15:26:11 (+0:06:40) ± 0:03:38"
    assert analysis.format_release_time(res, "%H:%M:%S") == "15:26:11 (+0:06:40) ± 0:03:38"
    assert analysis.format_error(0.123) == "0.12" and analysis.format_error(float("nan")) == "n/a"
    assert analysis.format_timedelta(pd.NaT) == "n/a"


# ---------------------------------------------------------------- reference times

def test_reference_time_at_sun_observed_at_earth():
    # Earth at 0.9926 AU on 2021-10-28: 8 min 15.8 s of light travel time
    reference = analysis.reference_time_at_sun(("SXR peak", datetime(2021, 10, 28, 15, 35)))
    assert reference.label == "SXR peak" and reference.observer == "Earth"
    assert reference.observed == pd.Timestamp("2021-10-28 15:35")
    assert (reference.observed - reference.at_sun).total_seconds() == pytest.approx(495.75, abs=0.1)
    assert analysis.reference_time_at_sun(("SXR peak", datetime(2021, 10, 28, 15, 35), "earth")) == reference


def test_reference_time_at_sun_of_a_distance_and_a_spacecraft(monkeypatch):
    reference = analysis.reference_time_at_sun(("type III", "2021-10-28 15:35", 0.5))
    assert reference.observer == "0.50 AU"
    assert (reference.observed - reference.at_sun).total_seconds() == pytest.approx(0.5 * 499.005, abs=1e-3)
    # the distance of a spacecraft at the observed time
    times = []
    monkeypatch.setattr(observers.StereoA, "distance", lambda self, time: times.append(time) or 0.96 * u.AU)
    reference = analysis.reference_time_at_sun(("type III", datetime(2021, 10, 28, 15, 35), "sta"))
    assert reference.observer == "STEREO-A" and times == [pd.Timestamp("2021-10-28 15:35")]
    assert (reference.observed - reference.at_sun).total_seconds() == pytest.approx(0.96 * 499.005, abs=1e-3)


@pytest.mark.parametrize("observer, message", [
    (0, r"Invalid distance of the observer: 0 AU. Give the distance from the Sun in AU \(e.g. 0.98\), above 0 and at most 10"),
    (-1.0, "Invalid distance of the observer: -1 AU"),
    (1.496e8, "Invalid distance of the observer: 1.496e[+]08 AU"),
    ("goes", r"Unknown observer 'goes': use 'earth', 'solo' \(Solar Orbiter\), 'sta' \(STEREO-A\), 'psp' \(Parker Solar "
             r"Probe\), 'soho' \(SOHO\), 'wind' \(Wind\), 'bepi' \(BepiColombo\), or a distance in AU"),
])
def test_reference_time_at_sun_invalid_observer(observer, message):
    with pytest.raises(ValueError, match=message):
        analysis.reference_time_at_sun(("SXR peak", datetime(2021, 10, 28, 15, 35), observer))


def test_vda_results_at_the_sun_at_1_au_and_at_the_spacecraft(fixed_spacecraft_distance):
    v = vda_with_onsets(["2021-10-31 03:20", "2021-10-31 03:40", "2021-10-31 04:00", "2021-10-31 04:21"])
    v.compute_vda()
    assert list(v.results.columns) == ["Release Time (Sun)", "Release Time (1 AU)", "Release Time", "Release Time Error",
                                       "Extra Time", "APL", "APL Error"]
    res = v.results.loc[1]
    at_sun, at_au = analysis.release_times(v.vda_fits[1])
    assert res["Release Time (Sun)"] == f"{at_sun:%Y-%m-%d %H:%M:%S}"
    assert res["Release Time (1 AU)"] == f"{at_au:%Y-%m-%d %H:%M:%S}"
    # the time at the spacecraft adds its light travel time from the Sun
    assert pd.Timestamp(res["Release Time"]) - pd.Timestamp(res["Release Time (Sun)"]) == pytest.approx(
        res["Extra Time"], abs=pd.Timedelta(seconds=1))


@pytest.mark.parametrize("observer", ["sta", "STEREO-A", "stereo-a", " Stereo-A "])
def test_reference_time_observer_by_name_or_label(monkeypatch, observer):
    monkeypatch.setattr(observers.StereoA, "distance", lambda self, time: 0.96 * u.AU)
    reference = analysis.reference_time_at_sun(("type III", datetime(2021, 10, 28, 15, 35), observer))
    assert reference.observer == "STEREO-A"
    assert (reference.observed - reference.at_sun).total_seconds() == pytest.approx(0.96 * 499.005, abs=1e-3)


def test_reference_time_observed_at_earth_by_label():
    assert analysis.reference_time_at_sun(("SXR", datetime(2021, 10, 28, 15, 35), "Earth")) == \
        analysis.reference_time_at_sun(("SXR", datetime(2021, 10, 28, 15, 35)))
