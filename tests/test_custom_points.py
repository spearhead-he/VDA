"""Custom points: the VDA fit of onset points given by the user, its points file and its notebook"""
from datetime import datetime

import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from ipywidgets import widgets

from spearhead.vda import notebook, views
from spearhead.vda.analysis import (CustomPointsVDA, inverse_beta, mean_energy, parse_energy, parse_particle,
                                    points_df, read_points_file)
from spearhead.vda.notebook import VDA_custom_points_notebook
from tests.helpers import write_csv
from tests.test_notebook import find

EXAMPLE = "examples/custom_points_example.csv"
# light travel time from the Sun at 1 AU, in seconds
AU_SECONDS = (1 * u.AU / const.c).to(u.s).value


@pytest.fixture
def displayed(monkeypatch):
    objects = []
    monkeypatch.setattr(notebook, "display", objects.append)
    yield objects
    plt.close("all")


def line_points(particle="protons"):
    """Points on a line: 1200 s per inverse beta, crossing inverse beta 0 at 03:00, with the energies of the
    inverse betas 1.5, 2, 3 and 4"""
    inv_betas = [1.5, 2, 3, 4]
    rest = 938.27 if particle == "protons" else 0.511
    energies = [rest * (1 / np.sqrt(1 - 1 / b ** 2) - 1) for b in inv_betas]
    onsets = [datetime(2021, 10, 31, 3) + pd.Timedelta(seconds=1200 * b) for b in inv_betas]
    return [(onset, f"{energy}", particle) for onset, energy in zip(onsets, energies)]


# ---------------------------------------------------------------- energies and particles

def test_inverse_beta():
    assert inverse_beta(mean_energy(10, 40), "protons") == pytest.approx(1 / np.sqrt(1 - (938.27 / 958.27) ** 2))
    assert mean_energy(10, 40) == pytest.approx(20)
    # electrons are faster than protons of the same energy
    assert inverse_beta(1, "electrons") < inverse_beta(1, "protons")


@pytest.mark.parametrize("text, energies", [("13.5", (13.5, 13.5)), (" 13-16 ", (13, 16)), (0.7, (0.7, 0.7))])
def test_parse_energy(text, energies):
    assert parse_energy(text) == energies


@pytest.mark.parametrize("text", ["", "a", "1-2-3", "16-13", "0", "-5"])
def test_parse_invalid_energy(text):
    with pytest.raises(ValueError, match="Invalid energy"):
        parse_energy(text)


def test_parse_particle():
    assert [parse_particle(p) for p in ("Protons", "p", " electrons", "E")] == [
        "protons", "protons", "electrons", "electrons"]
    with pytest.raises(ValueError, match="Unknown particle 'ions'"):
        parse_particle("ions")


def test_points_df():
    df = points_df([("2024-11-21 18:20:14", "0.67-10.04", "e"), (datetime(2024, 11, 21, 18, 28), "64", "protons")])
    assert list(df.columns) == ["Onset Time", "Particle", "Low Energy", "High Energy"]
    assert df["Particle"].tolist() == ["electrons", "protons"]
    assert df["Low Energy"].tolist() == [0.67, 64] and df["High Energy"].tolist() == [10.04, 64]
    with pytest.raises(ValueError, match="Point 2: Invalid energy"):
        points_df([("2024-11-21 18:20", "1", "p"), ("2024-11-21 18:30", "x", "p")])


# ---------------------------------------------------------------- VDA fit

def test_custom_points_exact_line():
    vda = CustomPointsVDA(points_df(line_points()), 1.0)
    vda.compute_vda()
    result = vda.results.loc[1]
    assert pd.Timestamp(result["Release Time"]) == pd.Timestamp("2021-10-31 03:00") + pd.Timedelta(seconds=round(AU_SECONDS))
    assert result["APL"] == pytest.approx(1200 / AU_SECONDS)
    assert pd.Timedelta(result["Extra Time"]).total_seconds() == pytest.approx(AU_SECONDS)
    fit = vda.vda_fits[1]
    assert fit["inv_betas"] == pytest.approx([1.5, 2, 3, 4])
    assert list(fit["particles"]) == ["protons"] * 4


def test_custom_points_sorted_by_inverse_beta():
    # protons of inverse beta 4, 3 and 2, and electrons of inverse beta 1.5
    vda = CustomPointsVDA(points_df(line_points()[:0:-1] + line_points("electrons")[:1]), 0.5)
    vda.compute_vda()
    fit = vda.vda_fits[1]
    assert list(fit["inv_betas"]) == sorted(fit["inv_betas"])
    assert list(fit["particles"]) == ["electrons", "protons", "protons", "protons"]


def test_custom_points_need_two_points(capsys):
    vda = CustomPointsVDA(points_df(line_points()[:1]), 1.0, name="group SOHO")
    vda.compute_vda()
    assert "Not enough onset points in group SOHO" in capsys.readouterr().out
    assert vda.results.loc[1].isna().all() and vda.vda_fits == {}


def test_custom_points_published_path_length():
    # the onsets of the 10 November 2022 event, Solar Orbiter EPT-South electrons at 0.61 AU, read from the VDA figure
    # of Rodriguez-Garcia et al. (arXiv:2604.12731): L = 0.66 +/- 0.08 AU, release at the Sun 17:08:43
    energies = [0.036, 0.041, 0.047, 0.054, 0.063, 0.074, 0.087, 0.102, 0.121, 0.143, 0.169, 0.192]
    onsets = ["17:24:00", "17:23:30", "17:23:00", "17:21:00", "17:19:30", "17:20:30", "17:19:00", "17:19:00",
              "17:18:00", "17:17:30", "17:18:00", "17:16:00"]
    vda = CustomPointsVDA(points_df([(f"2022-11-10 {onset}", f"{energy}", "e") for onset, energy in zip(onsets, energies)]),
                          0.61)
    vda.compute_vda()
    assert vda.results.loc[1, "APL"] == pytest.approx(0.66, abs=0.02)
    fit = vda.vda_fits[1]
    release_at_the_sun = fit["t0"] + pd.Timedelta(seconds=fit["b"])
    assert abs(release_at_the_sun - pd.Timestamp("2022-11-10 17:08:43")) < pd.Timedelta(seconds=30)


# ---------------------------------------------------------------- points file

def test_read_example_points_file():
    groups = read_points_file(EXAMPLE)
    assert list(groups) == ["SOHO", "STEREO-A", "Solar Orbiter"]
    assert [len(g.points) for g in groups.values()] == [10, 13, 11]
    assert groups["Solar Orbiter"].distance == 0.8269
    assert groups["SOHO"].name == "group SOHO"


def test_points_file_without_groups(tmp_path):
    path = write_csv(tmp_path, "points.csv", "Distance,Particle,Energy,Onset Time\n"
                                             "0.9,p,10,2024-11-21 18:30\n0.9,e,1-2,2024-11-21 18:20\n")
    groups = read_points_file(path)
    assert list(groups) == ["Points"]


@pytest.mark.parametrize("text, message", [
    ("Group,Particle,Energy,Onset Time\nA,p,10,2024-11-21 18:30\n", "Missing columns .*: Distance"),
    ("Group,Distance,Particle,Energy,Onset Time\nA,0.9,p,10,2024-11-21 18:30\nA,0.8,p,20,2024-11-21 18:20\n",
     "Group 'A' .*: give one distance"),
    ("Group,Distance,Particle,Energy,Onset Time\nA,0.9,p,10,2024-11-21 18:30\nA,0.9,x,20,2024-11-21 18:20\n",
     "Group 'A' .*: Point 2: Unknown particle 'x'"),
])
def test_invalid_points_file(tmp_path, text, message):
    with pytest.raises(ValueError, match=message):
        read_points_file(write_csv(tmp_path, "points.csv", text))


# ---------------------------------------------------------------- plot

def test_plot_of_custom_points(tmp_path):
    vdas = read_points_file(EXAMPLE)
    for vda in vdas.values():
        vda.compute_vda()
    fig = views.plot_vda_combined(vdas, 1, filename=tmp_path / "points.png", title="SEP event of 2024-11-21")
    assert fig.get_suptitle() == "SEP event of 2024-11-21"
    assert [line.get_label().split(":")[0] for line in fig.axes[0].get_lines()] == ["SOHO", "STEREO-A", "Solar Orbiter"]
    assert (tmp_path / "points.png").stat().st_size > 0
    assert views.custom_points_plot_filename(vdas) == "2024-11-21_1809_custom_points.png"
    plt.close(fig)


# ---------------------------------------------------------------- notebook

def test_parameters_form(displayed):
    tool = VDA_custom_points_notebook()
    tabs = tool.parameters_form()
    assert [tabs.get_title(i) for i in range(len(tabs.children))] == ["Points", "Views/Plots"]
    # a single group, which cannot be removed
    assert tool.groups == [{"label": "Group 1", "distance": 1.0, "points": []}]
    assert find(tabs, description="Remove Group")[0].disabled
    find(tabs.children[1], widgets.Text, "Title:")[0].value = "My event"
    assert tool.title == "My event"


def test_groups_widgets(displayed):
    tool = VDA_custom_points_notebook()
    tool.load_points(EXAMPLE)
    assert [g["label"] for g in tool.groups] == ["SOHO", "STEREO-A", "Solar Orbiter"]
    assert tool.groups[0]["points"][0] == ["2024-11-21 18:20:14", "0.67-10.04", "electrons"]
    assert tool.groups[2]["points"][0] == ["2024-11-21 18:34:34", "63.1-81.64", "protons"]
    tab = tool._points_tab()
    group_tabs = find(tab, widgets.Tab)[0]
    assert [group_tabs.get_title(i) for i in range(3)] == ["SOHO", "STEREO-A", "Solar Orbiter"]

    # the label is the title of the tab
    soho = group_tabs.children[0]
    find(soho, widgets.Text, "Label:")[0].value = "SOHO/ERNE"
    assert tool.groups[0]["label"] == "SOHO/ERNE" and group_tabs.get_title(0) == "SOHO/ERNE"
    find(soho, description="Distance from the Sun (AU):")[0].value = 0.98
    assert tool.groups[0]["distance"] == 0.98
    # the points: add (a copy of the last one), change and remove
    find(soho, description="Add Point")[0].click()
    assert len(tool.groups[0]["points"]) == 11
    assert tool.groups[0]["points"][10] == tool.groups[0]["points"][9]
    find(soho, widgets.Text, "Point 11 onset:")[0].value = "2024-11-21 19:10"
    find(soho, widgets.Text, "energy (MeV):")[10].value = "10-13"
    find(soho, widgets.Dropdown, "particle:")[10].value = "electrons"
    assert tool.groups[0]["points"][10] == ["2024-11-21 19:10", "10-13", "electrons"]
    find(soho, description="Remove")[0].click()
    assert len(tool.groups[0]["points"]) == 10
    assert tool.groups[0]["points"][0][0] == "2024-11-21 18:18:14"

    # groups: add and remove
    find(tab, description="Add Group")[0].click()
    assert tool.groups[3] == {"label": "Group 4", "distance": 1.0, "points": []}
    assert group_tabs.selected_index == 3
    find(group_tabs.children[3], description="Add Point")[0].click()
    assert tool.groups[3]["points"] == [["", "", "protons"]]
    find(group_tabs.children[1], description="Remove Group")[0].click()
    assert [g["label"] for g in tool.groups] == ["SOHO/ERNE", "Solar Orbiter", "Group 4"]
    assert [group_tabs.get_title(i) for i in range(3)] == ["SOHO/ERNE", "Solar Orbiter", "Group 4"]


def test_points_file_grays_out_the_groups(displayed, tmp_path):
    tool = VDA_custom_points_notebook()
    tab = tool._points_tab()
    wgt_file = find(tab, widgets.Text, "Points file:")[0]
    wgt_info = tab.children[0].children[2]
    groups = tab.children[1]

    def groups_disabled():
        return {w.disabled for w in find(groups) if hasattr(w, "disabled")}

    assert groups_disabled() == {False, True}
    wgt_file.value = str(tmp_path / "missing.csv")
    assert tool.points_filepath == str(tmp_path / "missing.csv")
    assert groups_disabled() == {True}
    assert "File not found" in wgt_info.value
    wgt_file.value = EXAMPLE
    assert wgt_info.value == "<i>3 groups, 34 points: SOHO, STEREO-A, Solar Orbiter</i>"
    wgt_file.value = write_csv(tmp_path, "bad.csv", "Particle,Energy\np,10\n")
    assert "Missing columns" in wgt_info.value
    wgt_file.value = ""
    assert wgt_info.value == ""
    # the remove button of a single group stays disabled
    assert find(groups, description="Remove Group")[0].disabled


def test_default_reference_time(displayed):
    tool = VDA_custom_points_notebook()
    tool._plots_tab()
    # without points, the current time
    assert abs(tool._default_reference_time() - datetime.now()).total_seconds() < 120
    tool.load_points(EXAMPLE)
    tab = tool._plots_tab()
    find(tab, description="Add Reference Time")[0].click()
    # the earliest onset, to the minute
    assert tool.reference_times == [("Reference 1", datetime(2024, 11, 21, 18, 9))]


def test_run_vda_from_the_form(displayed, monkeypatch, tmp_path, capsys):
    monkeypatch.chdir(tmp_path)
    tool = VDA_custom_points_notebook()
    tool.groups = [{"label": "A", "distance": 1.0, "points": line_points()},
                   {"label": "B", "distance": 0.5, "points": line_points("electrons")}]
    tool.reference_times = [("SXR peak", datetime(2021, 10, 31, 3, 5)), ("Not set", None)]
    tool.run_vda()
    out = capsys.readouterr().out
    assert "A\n    Release Time :" in out and "B\n    Release Time :" in out
    assert list(tool.results.index) == ["A", "B"] and tool.results.index.name == "Group"
    assert tool.results.loc["A", "APL"] == pytest.approx(1200 / AU_SECONDS)
    # the same slope at another distance: the same path length
    assert tool.results.loc["B", "APL"] == pytest.approx(1200 / AU_SECONDS)
    assert [type(o).__name__ for o in displayed] == ["Image"]
    assert [f.name for f in tmp_path.glob("*.png")] == ["2021-10-31_0330_custom_points.png"]


def test_run_vda_from_the_points_file(displayed, capsys):
    tool = VDA_custom_points_notebook()
    tool.save_plots = False
    tool.points_filepath = EXAMPLE
    tool.run_vda()
    assert list(tool.results.index) == ["SOHO", "STEREO-A", "Solar Orbiter"]
    assert tool.results.loc["STEREO-A", "Release Time"] == "2024-11-21 18:05:27"
    assert tool.results.loc["Solar Orbiter", "APL"] == pytest.approx(1.87, abs=0.005)


@pytest.mark.parametrize("groups, message", [
    ([{"label": "A", "distance": 1.0, "points": [("2021-10-31 03:30", "x", "p")]}], "Group 'A': Point 1: Invalid energy"),
    ([{"label": "A", "distance": 1.0, "points": []}, {"label": "A", "distance": 1.0, "points": []}],
     "Group 'A': give each group a different label"),
    ([{"label": " ", "distance": 1.0, "points": []}], "Group '': give each group a different label"),
])
def test_run_vda_invalid_groups(displayed, groups, message):
    tool = VDA_custom_points_notebook()
    tool.groups = groups
    with pytest.raises(ValueError, match=message):
        tool.run_vda()
