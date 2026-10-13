"""VDA fit and plot"""
import os
import time
from datetime import datetime, timedelta

import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
import pytest

from spearhead.vda import observers, views
from tests.helpers import make_vda


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
    points = {c.get_label(): len(c.get_offsets()) for c in fig.axes[0].collections}
    assert points["Protons"] == 2 and points["Electrons"] == 2
    plt.close(fig)
