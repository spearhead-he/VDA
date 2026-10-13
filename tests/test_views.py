"""Plots: background window, detected onsets and VDA fit, checked by their contents"""
from datetime import datetime

import astropy.units as u
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from spearhead.vda import observers, views
from tests.helpers import epd_like_data, make_vda

START, END = datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 20)
RISE = pd.Timestamp("2021-10-28 16:00")


def rising_epd_data(**kwargs):
    """Data as returned by epd_load, with a rise at 16:00 in the sun viewing only"""
    index = pd.date_range(START, END, freq="1min")
    df_protons, df_electrons, energies = epd_like_data(index, kwargs["sensor"])
    noise = np.where(np.arange(len(index)) % 2 == 0, 1.0, 1.2)
    rise = np.where(index >= RISE, 100.0, 0.0) if kwargs["viewing"] == "sun" else 0.0
    for df in (df_protons, df_electrons):
        for column in df.columns:
            df[column] = noise + rise
    return df_protons, df_electrons, energies


@pytest.fixture
def vda_onsets(monkeypatch):
    """Returns a function giving the VDA of the event with its onsets, from the synthetic data"""
    monkeypatch.setattr(observers, "epd_load", rising_epd_data)

    def onsets(viewings=("sun", "north"), onset_method="sigma"):
        v = make_vda(
            date_ranges=[(START, END)], resample_frequency="5min", viewings=list(viewings), onset_method=onset_method,
            channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}},
                            "electrons": {"HET/electrons Channel 1": {"sensor": "het", "channels": [0, 1]}}},
        )
        v.construct_times_df()
        v.construct_energies_df()
        v.construct_particles_df()
        v.group_energy_channels()
        v.calculate_onsets()
        v.clean_onsets()
        return v

    yield onsets
    plt.close("all")


def onset_figure(v, **kwargs):
    return views.plot_onsets(v, 1, "het", "protons", "H_Flux", "HET/protons Channel 1", **kwargs)


def as_datetime(x):
    """Datetime of an x value of a plot, given as a date number or a timestamp"""
    return pd.Timestamp(mdates.num2date(x) if isinstance(x, (float, np.floating)) else x).tz_localize(None)


def legend_labels(ax):
    legend = ax.get_legend()
    return [t.get_text() for t in legend.get_texts()] if legend else []


# ---------------------------------------------------------------- background window

def test_plot_bg(vda_onsets):
    v = vda_onsets()
    fig = views.plot_bg(v, 1)
    ax = fig.axes[0]
    # one line per grouped channel and viewing, with its energy range
    assert [line.get_label() for line in ax.get_lines()] == [
        "HET protons 2.00-8.00 MeV (sun)", "HET protons 2.00-8.00 MeV (north)",
        "HET electrons 0.50-2.00 MeV (sun)", "HET electrons 0.50-2.00 MeV (north)",
    ]
    assert [line.get_linestyle() for line in ax.get_lines()] == ["-", "--", "-", "--"]
    # protons and electrons with different colormaps
    colors = [mcolors.to_hex(line.get_color()) for line in ax.get_lines()]
    assert colors[0] != colors[2]
    assert ax.get_yscale() == "log"
    assert legend_labels(ax)[-1] == "Background"
    # the default background window: the first hour of the data
    span = [p for p in ax.patches if p.get_label() == "Background"][0]
    assert as_datetime(span.get_x()) == START
    assert as_datetime(span.get_x() + span.get_width()) == datetime(2021, 10, 28, 15)
    assert ax.get_title() == "Event 1\nBackground: 2021-10-28 14:00 to 15:00 (default)"


def test_plot_bg_set_window_on_another_day(vda_onsets):
    v = vda_onsets()
    v.set_bg_window(1, "2021-10-28 18:00", "2021-10-29 01:00")
    fig = views.plot_bg(v, 1)
    assert fig.axes[0].get_title() == "Event 1\nBackground: 2021-10-28 18:00 to 2021-10-29 01:00 (set)"


def test_plot_bg_without_viewing_labels(vda_onsets):
    # with a single viewing, the labels do not name it
    v = vda_onsets(viewings=["sun"])
    fig = views.plot_bg(v, 1)
    assert [line.get_label() for line in fig.axes[0].get_lines()] == [
        "HET protons 2.00-8.00 MeV", "HET electrons 0.50-2.00 MeV"]


# ---------------------------------------------------------------- detected onsets

def test_plot_onsets_sigma(vda_onsets):
    v = vda_onsets()
    fig = onset_figure(v, selected_viewing="sun")
    ax_sun, ax_north, ax_legend = fig.axes
    assert ax_sun.get_title() == "sun (selected)"
    assert ax_north.get_title() == "north"
    assert [line.get_label() for line in ax_sun.get_lines()] == ["Data", "BG level", "Threshold", "Onset"]
    # the onset is the last point before the rise
    onset_line = ax_sun.get_lines()[-1]
    assert as_datetime(onset_line.get_xdata()[0]) == RISE - pd.Timedelta("5min")
    method_specific = v.df_onsets_existing.loc[(1, "het", "protons", "sun", "H_Flux", "HET/protons Channel 1"),
                                               "Method Specific"]
    assert ax_sun.get_lines()[2].get_ydata()[0] == pytest.approx(method_specific["threshold"])
    values = ax_sun.texts[0].get_text().split("\n")
    assert values[0] == "Onset 15:55"
    assert values[1].startswith("BG ") and values[2].startswith("Threshold ")
    # the north viewing has no rise
    assert [t.get_text() for t in ax_north.texts] == ["No onset"]
    assert not ax_north.get_lines()
    # one legend for all the viewings
    assert legend_labels(ax_legend) == ["Data", "BG sample", "BG level", "Threshold", "Onset"]
    assert fig.get_suptitle() == "Detected onsets for event 1 (2021-10-28) | het/protons (2.00-8.00 MeV)"


def test_plot_onsets_poisson_cusum(vda_onsets):
    v = vda_onsets(onset_method="poisson_cusum")
    fig = onset_figure(v)
    ax_sun = fig.axes[0]
    assert ax_sun.get_title() == "sun"
    assert [line.get_label() for line in ax_sun.get_lines()] == ["Data", "BG level", "μd", "Onset"]
    values = ax_sun.texts[0].get_text().split("\n")
    assert values[2].startswith("μd ") and values[3].startswith("k ") and ", h " in values[3]


def test_plot_onsets_without_any_onset(vda_onsets):
    v = vda_onsets(viewings=["north"])
    fig = onset_figure(v)
    ax_north, ax_legend = fig.axes
    assert [t.get_text() for t in ax_north.texts] == ["No onset"]
    assert ax_legend.get_legend() is None


def test_plot_onsets_grid(vda_onsets):
    # five viewings: two rows of three subplots, the last one empty, and the row of the legend
    v = vda_onsets(viewings=["sun", "asun", "north", "south", "omni"])
    fig = onset_figure(v)
    assert len(fig.axes) == 7
    assert [ax.get_title() for ax in fig.axes[:5]] == ["sun", "asun", "north", "south", "omni"]
    assert not fig.axes[5].axison and not fig.axes[6].axison


# ---------------------------------------------------------------- files

def test_plots_saved_to_files(vda_onsets, monkeypatch, tmp_path):
    v = vda_onsets()
    views.plot_bg(v, 1, filename=tmp_path / "bg.png")
    onset_figure(v, filename=tmp_path / "onsets.png")
    monkeypatch.setattr(observers.SolarOrbiter, "distance", lambda self, time: 0.8 * u.AU)
    v.construct_options_df()
    v.select_onsets()
    v.construct_energy_channels_characteristics()
    v.compute_vda()
    views.plot_vda(v, 1, filename=tmp_path / "vda.png")
    for name in ("bg.png", "onsets.png", "vda.png"):
        assert (tmp_path / name).stat().st_size > 0


@pytest.mark.parametrize("resample_frequency, suffix", [("5min", "5min"), ("", "noresample")])
def test_vda_plot_filename(vda_onsets, resample_frequency, suffix):
    v = vda_onsets()
    v.parameters.resample_frequency = resample_frequency
    assert views.vda_plot_filename(v, 1) == f"2021-10-28_1400_2021-10-28_2000_het-protons_het-electrons_{suffix}.png"
