"""Multi-spacecraft notebook: the parameters form, the analysis of each spacecraft and the combined VDA plot"""
from datetime import datetime

import astropy.units as u
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from ipywidgets import widgets

from spearhead.vda import OnsetSelection, VDA, analysis, notebook, observers, views
from spearhead.vda.notebook import VDA_multi_spacecraft_notebook
from tests.helpers import legend_texts, reference_legend, rising_epd_data, rising_observer_load, table_rows
from tests.test_notebook import find

START, END = datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 20)


@pytest.fixture
def displayed(monkeypatch):
    """The objects displayed by the notebook, with synthetic data of Solar Orbiter and STEREO-A, and no SOHO data"""
    def no_data(*args, **kwargs):
        raise ValueError("No SOHO ERNE data")

    monkeypatch.setattr(observers, "epd_load", rising_epd_data)
    monkeypatch.setattr(observers.StereoA, "load", rising_observer_load)
    monkeypatch.setattr(observers.Soho, "load", no_data)
    monkeypatch.setattr(VDA, "define_spacecraft_parameters", lambda self: None)
    monkeypatch.setattr(observers.SolarOrbiter, "distance", lambda self, time: 0.8 * u.AU)
    monkeypatch.setattr(observers.StereoA, "distance", lambda self, time: 0.96 * u.AU)
    objects = []
    monkeypatch.setattr(notebook, "display", objects.append)
    yield objects
    plt.close("all")


def make_tool(spacecraft=("solo", "sta")):
    """Notebook of the default event, with grouped channels of the synthetic data"""
    tool = VDA_multi_spacecraft_notebook(spacecraft)
    tool.parameters.view_dfs = False
    tool.save_plots = False
    tool.notebook("solo").parameters.channel_groups = {
        "protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}},
        "electrons": {"HET/electrons Channel 1": {"sensor": "het", "channels": [0, 1]}},
    }
    tool.notebook("sta").parameters.channel_groups = {
        "protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [0, 1]},
                    "HET/protons Channel 2": {"sensor": "het", "channels": [4, 5]}},
        "electrons": {"HET/electrons Channel 1": {"sensor": "het", "channels": [0, 1]}},
    }
    return tool


# ---------------------------------------------------------------- parameters form

def test_parameters_form(displayed):
    tool = make_tool()
    tabs = tool.parameters_form()
    assert [tabs.get_title(i) for i in range(len(tabs.children))] == ["Event", "Spacecraft", "Onsets", "Views/Plots"]
    event_tab, spacecraft_tab = tabs.children[:2]
    # the datetime range of the event
    find(event_tab, description="end:")[0].value = datetime(2021, 10, 28, 21)
    assert tool.parameters.date_ranges == [(START, datetime(2021, 10, 28, 21))]
    # the folder of the downloaded data
    find(event_tab, widgets.Text, "Data folder:")[0].value = "sep_data"
    assert tool.parameters.data_path == "sep_data"
    # a tab with the settings of each selected spacecraft
    settings = find(spacecraft_tab, widgets.Tab)[0]
    assert [settings.get_title(i) for i in range(len(settings.children))] == ["Solar Orbiter", "STEREO-A"]
    find(spacecraft_tab, widgets.Checkbox, "Solar Orbiter")[0].value = False
    assert tool.spacecraft == ["sta"]
    assert [settings.get_title(i) for i in range(len(settings.children))] == ["STEREO-A"]
    find(spacecraft_tab, widgets.Checkbox, "Solar Orbiter")[0].value = True
    # the spacecraft keep the order of OBSERVERS, and their settings
    assert tool.spacecraft == ["solo", "sta"]
    assert settings.get_title(0) == "Solar Orbiter"
    # the resample frequency of each spacecraft
    find(settings.children[1], widgets.Text, "Resample frequency:")[0].value = "1min"
    assert tool.notebook("sta").parameters.resample_frequency == "1min"
    assert tool.notebook("solo").parameters.resample_frequency == "5min"
    # the plot of each spacecraft is optional
    plots_tab = tabs.children[3]
    find(plots_tab, description="Show the VDA plot of each spacecraft, besides the combined plot")[0].value = True
    assert tool.show_spacecraft_plots is True


def test_select_unknown_spacecraft(displayed):
    with pytest.raises(ValueError, match="Unknown observer 'juice'"):
        VDA_multi_spacecraft_notebook(["juice"])


def test_shared_parameters(displayed, tmp_path):
    tool = make_tool()
    tool.parameters.date_ranges = [(START, datetime(2021, 10, 28, 19))]
    tool.parameters.onset_method = "poisson_cusum"
    tool.parameters.onset_method_parameters["cusum_minutes"] = 20
    tool.vda.set_default_bg_window(10, 50)
    tool.parameters.data_path = str(tmp_path)
    tool.run_data()
    for name in ("solo", "sta"):
        parameters = tool.notebook(name).parameters
        assert parameters.date_ranges == [(START, datetime(2021, 10, 28, 19))]
        assert parameters.data_path == str(tmp_path)
        assert parameters.onset_method == "poisson_cusum"
        assert parameters.onset_method_parameters["cusum_minutes"] == 20
        assert parameters.bg_after_start == (10, 50)
    # the parameters are copies
    tool.parameters.onset_method_parameters["cusum_minutes"] = 40
    assert tool.notebook("sta").parameters.onset_method_parameters["cusum_minutes"] == 20


# ---------------------------------------------------------------- analysis steps

def test_spacecraft_without_data_left_out(displayed, capsys):
    tool = make_tool(("solo", "sta", "soho"))
    tool.run_data()
    assert list(tool.analyses) == ["solo", "sta"]
    assert "Left out, no data: No SOHO ERNE data" in capsys.readouterr().out


class FakeRequestsConnectionError(OSError):
    """Like requests.exceptions.ConnectionError, which is not the built-in ConnectionError"""


FakeRequestsConnectionError.__name__ = "ConnectionError"


@pytest.mark.parametrize("error", [ConnectionError("Connection aborted."), FakeRequestsConnectionError("Max retries"),
                                   TimeoutError("timed out")])
def test_spacecraft_without_connection_left_out(displayed, monkeypatch, capsys, error):
    def no_connection(*args, **kwargs):
        raise error

    monkeypatch.setattr(observers.StereoA, "load", no_connection)
    tool = make_tool()
    tool.run_data()
    assert list(tool.analyses) == ["solo"]
    assert ("Left out, the data could not be downloaded (the data server may be unavailable, run the data step "
            f"again later): {error}") in capsys.readouterr().out


def test_without_spacecraft(displayed):
    tool = make_tool()
    tool.select_spacecraft("solo", False)
    tool.select_spacecraft("sta", False)
    with pytest.raises(ValueError, match="No spacecraft"):
        tool.run_data()
    tool = make_tool(("soho",))
    with pytest.raises(ValueError, match="No spacecraft with data"):
        tool.run_data()


def test_background_and_onset_selection_per_spacecraft(displayed):
    tool = make_tool()
    tool.parameters.onset_selection = OnsetSelection.INTERACTIVE
    tool.run_data()
    box = tool.background_selection()
    wgt_spacecraft = box.children[0]
    assert wgt_spacecraft.options == (("Solar Orbiter", "solo"), ("STEREO-A", "sta"))
    wgt_spacecraft.value = "sta"
    # the background window of STEREO-A
    wgt_bg = find(box, widgets.SelectionRangeSlider)[0]
    wgt_bg.index = (12, 24)
    assert tool.notebook("sta").vda.bg_window_source(1) == "set"
    assert tool.notebook("solo").vda.bg_window_source(1) == "default"

    box = tool.run_onsets()
    box.children[0].value = "sta"
    wgt_channel = find(box, widgets.Dropdown, "Channel:")[0]
    assert [label for label, _ in wgt_channel.options] == [
        "HET protons 1.00-3.00 MeV", "HET protons 5.00-7.00 MeV", "HET electrons 0.10-0.30 MeV"]
    find(box, widgets.Dropdown, "Selected viewing:")[0].value = None
    assert tool.notebook("sta").parameters.selected_onsets["Viewing"].tolist() == [None, "omni", "omni"]


def test_spacecraft_without_onsets(displayed, capsys):
    tool = make_tool(("solo",))
    tool.notebook("solo").parameters.viewings = ["north"]
    tool.parameters.onset_selection = OnsetSelection.INTERACTIVE
    tool.parameters.view_dfs = True
    tool.run_data()
    capsys.readouterr()
    box = tool.run_onsets()
    # the name of the spacecraft before its onsets
    assert capsys.readouterr().out.startswith("Solar Orbiter:\nEvent 1: no onsets in its 2 grouped channels.")
    assert box.children[1].children[0].value == "No onsets were determined"
    tool.run_vda()
    assert "Solar Orbiter\n    No results (not enough onset points" in capsys.readouterr().out
    assert pd.isna(tool.results.loc["Solar Orbiter"]).all()


def test_run_vda(displayed, monkeypatch, tmp_path, capsys):
    monkeypatch.chdir(tmp_path)
    tool = make_tool()
    tool.save_plots = True
    tool.reference_times = [("SXR peak", datetime(2021, 10, 28, 15, 50)), ("Outside", datetime(2021, 10, 29))]
    tool.run_data()
    assert tool.run_onsets() is None
    tool.run_vda()
    out = capsys.readouterr().out
    assert "Solar Orbiter\n    Release at the Sun : " in out and "STEREO-A\n    Release at the Sun : " in out
    assert "\n    At the spacecraft  : " in out
    assert list(tool.results.index) == ["Solar Orbiter", "STEREO-A"]
    assert tool.results.index.name == "Spacecraft"
    assert list(tool.results.columns) == analysis.RESULTS_COLUMNS
    # only the combined plot, saved
    assert [type(o).__name__ for o in displayed] == ["Image"]
    assert [f.name for f in tmp_path.glob("*.png")] == ["2021-10-28_1400_2021-10-28_2000_solo-sta.png"]

    displayed.clear()
    tool.show_spacecraft_plots = True
    tool.parameters.view_dfs = True
    tool.run_vda()
    # the plot of each spacecraft, the combined plot and the results table
    assert [type(o).__name__ for o in displayed] == ["Image", "Image", "Image", "DataFrame"]
    assert sorted(f.name for f in tmp_path.glob("*.png")) == [
        "2021-10-28_1400_2021-10-28_2000_solo-sta.png",
        "solo_2021-10-28_1400_2021-10-28_2000_het-protons_het-electrons_5min.png",
        "sta_2021-10-28_1400_2021-10-28_2000_het-protons_het-electrons_5min.png",
    ]


# ---------------------------------------------------------------- combined plot

@pytest.fixture
def vdas(displayed):
    tool = make_tool()
    tool.run_data()
    tool.run_onsets()
    tool.run_vda()
    return {"Solar Orbiter": tool.notebook("solo").vda, "STEREO-A": tool.notebook("sta").vda}


def test_plot_vda_combined(vdas):
    references = [("SXR peak", datetime(2021, 10, 28, 15, 50))]
    fig = views.plot_vda_combined(vdas, 1, references=references)
    ax = fig.axes[0]
    assert fig.get_suptitle() == "VDA of 2021-10-28"
    # one fit line per spacecraft, in its color, from inverse beta 0
    lines = ax.get_lines()
    assert len(lines) == 2
    assert all(line.get_xdata()[0] == 0 for line in lines)
    assert mcolors.to_hex(lines[0].get_color()) != mcolors.to_hex(lines[1].get_color())
    # the symbols and reference times in the legend, the results of each spacecraft in the table, in its color
    assert legend_texts(ax) == ["Protons", "Electrons", "VDA fits", "2σ bands", "Release at the Sun"]
    references = reference_legend(ax)
    assert references.get_title().get_text() == ""
    assert [t.get_text() for t in references.get_texts()] == ["SXR peak: 15:41:44 (15:50:00 at Earth)"]
    # the legends inside the plot are a bit smaller than the legend font size
    assert ax.get_legend().get_texts()[0].get_fontsize() == pytest.approx(0.9 * plt.rcParams["legend.fontsize"])
    assert references.get_texts()[0].get_fontsize() == ax.get_legend().get_texts()[0].get_fontsize()
    rows = table_rows(fig)
    assert rows[0] == ["", "Release at the Sun", "At 1 AU", "Error", "APL (AU)"]
    assert [row[0] for row in rows[1:]] == ["Solar Orbiter", "STEREO-A"]
    label_colors = [mcolors.to_hex(fig.axes[1].tables[0][row, 0].get_text().get_color()) for row in (1, 2)]
    assert label_colors == [mcolors.to_hex(line.get_color()) for line in lines]
    assert ax.get_xlim()[0] == 0
    # the points of each spacecraft in its color
    points = [c for c in ax.collections if len(c.get_offsets()) and c.get_label() != "SXR peak"
              and c.get_offsets()[0][0] != 0]
    assert {mcolors.to_hex(c.get_facecolor()[0]) for c in points} == {mcolors.to_hex(line.get_color()) for line in lines}


def test_plot_vda_combined_without_fit(vdas):
    # a spacecraft without a fit of the event is left out
    vdas["SOHO"] = VDA(vdas["STEREO-A"].parameters)
    fig = views.plot_vda_combined(vdas, 1)
    assert [row[0] for row in table_rows(fig)[1:]] == ["Solar Orbiter", "STEREO-A"]
    fig = views.plot_vda_combined({"SOHO": vdas["SOHO"]}, 1)
    assert fig.axes[0].get_lines() == [] and fig.get_suptitle() == "" and len(fig.axes[1].tables) == 0


def test_vda_combined_plot_filename(vdas, tmp_path):
    assert views.vda_combined_plot_filename(vdas) == "2021-10-28_1400_2021-10-28_2000_solo-sta.png"
    views.plot_vda_combined(vdas, 1, filename=tmp_path / "combined.png")
    assert (tmp_path / "combined.png").stat().st_size > 0
    assert pd.notna(vdas["STEREO-A"].results.loc[1, "APL"])


def test_plot_vda_combined_release_times(vdas):
    ax = views.plot_vda_combined(vdas, 1).axes[0]
    assert ax.get_ylabel() == "Time (UTC)"
    # the release time at the Sun, where the fit line meets inverse beta 0, and at 1 AU
    fig = ax.figure
    at_sun, at_au = analysis.release_times(vdas["Solar Orbiter"].vda_fits[1])
    assert table_rows(fig)[1][1:3] == [f"{at_sun:%H:%M:%S}", f"{at_au:%H:%M:%S}"]
