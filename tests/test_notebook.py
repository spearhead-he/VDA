"""Notebook interface: the parameters form, the background window and onset selection, and the analysis steps"""
from datetime import datetime

import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from ipywidgets import widgets

from spearhead.vda import OnsetSelection, VDA, VDA_parameters
from spearhead.vda import notebook, observers
from spearhead.vda.notebook import VDA_notebook
from tests.helpers import epd_like_data, make_vda, rising_epd_data, write_csv


def test_form_without_data_of_the_observer(monkeypatch):
    # the error is shown in the form, which keeps working
    def no_data(*args, **kwargs):
        raise ValueError("No BepiColombo SIXS-P data")

    monkeypatch.setattr(observers.BepiColombo, "load", no_data)
    monkeypatch.setattr(observers, "epd_load", lambda **kwargs: epd_like_data(
        pd.date_range("2021-10-28 14:00", periods=4, freq="1min"), kwargs["sensor"]))
    tool = VDA_notebook()
    tool._data_tab()
    tool._channels_tab()
    assert tool._observer_error.value == ""
    tool.parameters.observer = "bepi"
    section = tool._channels_section()
    assert "color: red" in tool._observer_error.value and "No BepiColombo SIXS-P data" in tool._observer_error.value
    assert "color: red" in section.children[1].children[0].value
    assert tool.vda.df_energies is None
    # the grouped channels are grayed out
    wgt_channel_groups = section.children[-1]

    def disabled(widget):
        own = [widget.disabled] if hasattr(widget, "disabled") else []
        return own + [d for child in getattr(widget, "children", ()) for d in disabled(child)]

    assert disabled(wgt_channel_groups) and all(disabled(wgt_channel_groups))
    tool.parameters.observer = "solo"
    section = tool._channels_section()
    assert tool._observer_error.value == ""
    assert tool.vda.df_energies is not None
    assert not any(disabled(section.children[-1]))
    # the fake data have 3 proton channels: the other channels are listed without energy range
    assert tool._channel_options("het", "protons")[:4] == [("0: 1-2 MeV", 0), ("1: 2-4 MeV", 1), ("2: 4-8 MeV", 2), ("3", 3)]


def test_form_reads_the_channels_again(monkeypatch):
    # after an error, the button reads the channels again, e.g. after the events are changed
    calls = []

    def load(self, sensor, startdate, enddate, viewing, path, particles=None):
        calls.append(startdate)
        if len(calls) == 1:
            raise ValueError("No data")
        index = pd.date_range(startdate, periods=4, freq="1min")
        return {particle: observers.ParticleData(pd.DataFrame(1.0, index=index, columns=[1, 2]),
                                                 pd.DataFrame({"Low Energy": [1.0, 2.0], "Bin Width": [1.0, 1.0]},
                                                              index=[1, 2]))
                for particle in particles}

    monkeypatch.setattr(observers.BepiColombo, "load", load)
    tool = VDA_notebook(VDA(VDA_parameters(observer="bepi")))
    tool._data_tab()
    tool._channels_tab()
    assert "No data" in tool._observer_error.value
    error_row = tool._wrapper_channels.children[0].children[1]
    error_row.children[1].click()
    assert tool._observer_error.value == ""
    assert len(tool._wrapper_channels.children[0].children) == 2


# ---------------------------------------------------------------- parameters form

def find(widget, cls=None, description=None):
    """The widgets inside the widget (itself included) of the class and with the description"""
    found = []
    if cls is None or isinstance(widget, cls):
        if description is None or getattr(widget, "description", None) == description:
            found.append(widget)
    for child in getattr(widget, "children", ()):
        found += find(child, cls, description)
    return found


@pytest.fixture
def displayed(monkeypatch):
    """The objects displayed by the notebook (figures and the buttons of Add Channel), with synthetic data"""
    monkeypatch.setattr(observers, "epd_load", rising_epd_data)
    objects = []
    monkeypatch.setattr(notebook, "display", objects.append)
    yield objects
    plt.close("all")


def test_parameters_form_tabs(displayed):
    tabs = VDA_notebook().parameters_form()
    assert [tabs.get_title(i) for i in range(len(tabs.children))] == [
        "Events", "Data", "Energy channels", "Onsets", "Views/Plots"]


def test_events_file_grays_out_the_date_ranges(displayed, tmp_path):
    tool = VDA_notebook()
    tab = tool._events_tab()
    wgt_file = find(tab, widgets.Text)[0]
    wgt_ranges = tab.children[1]
    wrapper_file = tab.children[0].children[2]

    def ranges_disabled():
        return {w.disabled for w in find(wgt_ranges) if hasattr(w, "disabled")}

    # a single event cannot be removed
    assert find(wgt_ranges, description="Remove Event")[0].disabled
    assert not find(wgt_ranges, description="Add Event")[0].disabled

    wgt_file.value = str(tmp_path / "missing.csv")
    assert tool.parameters.input_filepath == str(tmp_path / "missing.csv")
    assert ranges_disabled() == {True}
    assert "File not found" in wrapper_file.children[0].value

    path = write_csv(tmp_path, "events.csv", "Event No,Start Time,BG End,End Time\n"
                                             "1,2021-05-22 19:45,2021-05-22 20:15,2021-05-23 02:45\n"
                                             "2,2021-11-09 15:25,,2021-11-09 22:25\n")
    wgt_file.value = path
    assert "2 events, 1 with a background window" in wrapper_file.children[0].value

    wgt_file.value = write_csv(tmp_path, "unknown.csv", "Event No,BG Start,End Time\n1,2021-05-22 19:45,2021-05-23 02:45\n")
    assert "Accepted columns" in wrapper_file.children[0].value

    # reference times: the hours prior to and after them
    wgt_file.value = write_csv(tmp_path, "reference.csv", "Event,Reference Time\n1,2021-05-22 21:45\n")
    sliders = find(wrapper_file, widgets.IntSlider)
    assert [s.description for s in sliders] == ["Hours prior to the reference time:", "Hours after the reference time:"]
    sliders[0].value = 3
    sliders[1].value = 7
    assert (tool.parameters.bg_hours_prior, tool.parameters.bg_hours_after) == (3, 7)

    wgt_file.value = ""
    assert wrapper_file.children == ()
    assert ranges_disabled() == {False, True}
    assert find(wgt_ranges, description="Remove Event")[0].disabled


def test_add_and_remove_events(displayed):
    tool = VDA_notebook()
    wgt_ranges = tool._date_ranges_widget()
    first = tool.parameters.date_ranges[0]
    find(wgt_ranges, description="Add Event")[0].click()
    # the new event is a copy of the last one
    assert tool.parameters.date_ranges == [first, first]
    pickers = find(wgt_ranges, widgets.widget_datetime.NaiveDatetimePicker)
    assert [p.description for p in pickers] == ["Event 1 start:", "end:", "Event 2 start:", "end:"]
    assert not any(b.disabled for b in find(wgt_ranges, description="Remove Event"))
    pickers[2].value = datetime(2022, 1, 1, 10)
    pickers[3].value = datetime(2022, 1, 1, 16)
    assert tool.parameters.date_ranges[1] == (datetime(2022, 1, 1, 10), datetime(2022, 1, 1, 16))
    find(wgt_ranges, description="Remove Event")[0].click()
    assert tool.parameters.date_ranges == [(datetime(2022, 1, 1, 10), datetime(2022, 1, 1, 16))]
    assert find(wgt_ranges, description="Remove Event")[0].disabled
    assert find(wgt_ranges, widgets.widget_datetime.NaiveDatetimePicker)[0].description == "Event 1 start:"


def test_data_tab_bindings(displayed):
    tool = VDA_notebook()
    tab = tool._data_tab()
    texts = {w.description: w for w in find(tab, widgets.Text)}
    texts["Resample frequency:"].value = "10min"
    texts["Load the data from:"].value = "saved.pkl"
    texts["Save the data to:"].value = "new.pkl"
    p = tool.parameters
    assert (p.resample_frequency, p.load_data_filepath, p.save_data_filepath) == ("10min", "saved.pkl", "new.pkl")


def test_viewings_checkboxes_keep_their_order(displayed):
    tool = VDA_notebook()
    tab = tool._data_tab()
    checkboxes = {w.description: w for w in find(tab, widgets.Checkbox)}
    assert list(checkboxes) == ["sun", "asun", "north", "south", "omni"]
    assert [w.value for w in checkboxes.values()] == [True, False, False, False, False]
    checkboxes["omni"].value = True
    checkboxes["north"].value = True
    assert tool.parameters.viewings == ["sun", "north", "omni"]
    checkboxes["sun"].value = False
    assert tool.parameters.viewings == ["north", "omni"]


def test_observer_change(displayed, monkeypatch):
    # the viewings and grouped channels are those of the new observer, with their notes
    def load(self, sensor, startdate, enddate, viewing, path, particles=None):
        index = pd.date_range(startdate, periods=4, freq="1min")
        channels = list(self.CHANNELS[sensor][particles[0]])
        return {particles[0]: observers.ParticleData(
            pd.DataFrame(1.0, index=index, columns=channels),
            pd.DataFrame({"Low Energy": [1.0] * len(channels), "Bin Width": [1.0] * len(channels)}, index=channels))}

    monkeypatch.setattr(observers.StereoA, "load", load)
    tool = VDA_notebook()
    tab = tool._data_tab()
    tool._channels_tab()
    wgt_observer = find(tab, widgets.Dropdown, "Spacecraft:")[0]
    assert wgt_observer.value == "solo"
    wgt_observer.value = "sta"
    assert tool.parameters.observer == "sta"
    viewings = tool._wrapper_viewings.children[0]
    assert [w.description for w in find(viewings, widgets.Checkbox)] == ["sun", "asun", "north", "south", "omni"]
    assert [w.description for w in find(viewings, widgets.Checkbox) if w.value] == ["sun", "omni"]
    assert "Viewings of each sensor: HET: omni; SEPT: sun, asun, north, south" in viewings.children[1].value
    labels = [w.value for w in find(tool._wrapper_channels, widgets.HTML) if "Channel" in w.value]
    assert len(labels) == len(tool.parameters.channel_groups["protons"]) + len(tool.parameters.channel_groups["electrons"])
    assert "STEREO-A" not in tool._observer_error.value


def test_add_and_remove_channel_groups(displayed):
    tool = VDA_notebook()
    tool.parameters.channel_groups = {"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}}}
    tool._data_tab()
    section = tool._channels_section()
    find(section, description="Add Channel")[0].click()
    # the buttons of the sensors and particles
    buttons = {b.description: b for b in find(displayed[-1], widgets.Button)}
    assert list(buttons) == ["HET/protons", "HET/electrons", "EPT/protons", "EPT/electrons"]
    buttons["HET/electrons"].click()
    buttons["HET/protons"].click()
    groups = tool.parameters.channel_groups
    assert list(groups["protons"]) == ["HET/protons Channel 1", "HET/protons Channel 2"]
    assert groups["electrons"] == {"HET/electrons Channel 1": {"sensor": "het", "channels": []}}

    selects = find(section, widgets.SelectMultiple)
    assert len(selects) == 3
    selects[1].value = (0, 1)
    assert groups["electrons"]["HET/electrons Channel 1"]["channels"] == [0, 1]

    # removing the last group of a particle removes the particle
    find(section, description="Remove Channel")[1].click()
    assert "electrons" not in groups
    assert len(find(section, widgets.SelectMultiple)) == 2
    find(section, description="Remove Channel")[0].click()
    assert list(groups["protons"]) == ["HET/protons Channel 2"]


def test_onset_method_change_resets_its_parameters(displayed):
    tool = VDA_notebook()
    tab = tool._onsets_tab()
    sliders = {w.description.split(" |")[0]: w for w in find(tab, widgets.IntSlider)}
    assert list(sliders) == ["s", "n"]
    sliders["s"].value = 4
    assert tool.parameters.onset_method_parameters == {"s": 4, "n": 3}
    find(tab, widgets.Dropdown, "Onset determination method:")[0].value = "poisson_cusum"
    assert tool.parameters.onset_method == "poisson_cusum"
    sliders = {w.description.split(" |")[0]: w for w in find(tab, widgets.IntSlider)}
    assert list(sliders) == ["cusum_minutes", "sigma_multiplier"]
    sliders["cusum_minutes"].value = 60
    assert tool.parameters.onset_method_parameters == {"cusum_minutes": 60, "sigma_multiplier": 2}


@pytest.mark.parametrize("resample_frequency, step", [("5min", 5), ("30s", 1), ("", 1)])
def test_default_background_slider(displayed, resample_frequency, step):
    tool = VDA_notebook(make_vda(resample_frequency=resample_frequency))
    slider = tool._bg_defaults_widget()
    assert slider.step == step
    # before the events are created, up to one day
    assert slider.max == 24 * 60
    slider.value = (10, 40)
    assert tool.parameters.bg_after_start == (10, 40)
    tool.vda.construct_times_df()
    assert tool._bg_defaults_widget().max == 6 * 60


def test_onset_selection_method_binding(displayed):
    tool = VDA_notebook()
    tab = tool._onsets_tab()
    find(tab, widgets.Dropdown, "Onset selection method:")[0].value = OnsetSelection.INTERACTIVE
    assert tool.parameters.onset_selection is OnsetSelection.INTERACTIVE


def test_plots_tab(displayed):
    tool = VDA_notebook()
    tab = tool._plots_tab()
    find(tab, description="Font size:")[0].value = 12
    find(tab, description="Legend font size:")[0].value = 10
    assert (plt.rcParams["font.size"], plt.rcParams["legend.fontsize"]) == (12, 10)
    find(tab, description="Save the VDA plots as .png files")[0].value = False
    assert tool.save_plots is False
    find(tab, description="Display the produced DataFrames")[0].value = True
    assert tool.parameters.view_dfs is True


# ---------------------------------------------------------------- analysis steps

def two_events_tool(**parameters):
    """Notebook of two events of the synthetic data, with one grouped channel of protons and one of electrons"""
    return VDA_notebook(make_vda(
        date_ranges=[(datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 20)),
                     (datetime(2021, 11, 9, 15), datetime(2021, 11, 9, 21))],
        resample_frequency="5min", viewings=["sun", "north"],
        channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}},
                        "electrons": {"HET/electrons Channel 1": {"sensor": "het", "channels": [0, 1]}}},
        **parameters,
    ))


def test_background_selection(displayed):
    tool = two_events_tool()
    tool.run_data()
    box = tool.background_selection()
    wgt_event = find(box, widgets.Dropdown, "Event:")[0]
    wgt_bg = find(box, widgets.SelectionRangeSlider)[0]
    wgt_status = box.children[2]
    assert [label for label, _ in wgt_event.options] == ["Event 1 (2021-10-28 14:00)", "Event 2 (2021-11-09 15:00)"]
    # the default window: the first hour, 13 points at 5 minutes
    assert wgt_bg.value == (pd.Timestamp("2021-10-28 14:00"), pd.Timestamp("2021-10-28 15:00"))
    assert wgt_status.value.startswith("Background 2021-10-28 14:00 to 2021-10-28 15:00 (default), 13 points")
    n_figures = len(displayed)

    wgt_bg.index = (12, 24)
    assert tool.vda.bg_window(1) == (pd.Timestamp("2021-10-28 15:00"), pd.Timestamp("2021-10-28 16:00"))
    assert tool.vda.bg_window_source(1) == "set"
    assert len(displayed) == n_figures + 1

    wgt_bg.index = (24, 24)
    assert wgt_status.value == "<b>The background start must be before its end</b>"
    assert tool.vda.bg_window(1)[1] == pd.Timestamp("2021-10-28 16:00")

    wgt_event.value = 2
    assert wgt_bg.value == (pd.Timestamp("2021-11-09 15:00"), pd.Timestamp("2021-11-09 16:00"))
    wgt_event.value = 1
    assert wgt_bg.value == (pd.Timestamp("2021-10-28 15:00"), pd.Timestamp("2021-10-28 16:00"))
    find(box, description="Reset to default")[0].click()
    assert tool.vda.bg_window_source(1) == "default"
    assert wgt_bg.value == (pd.Timestamp("2021-10-28 14:00"), pd.Timestamp("2021-10-28 15:00"))


def test_interactive_onset_selection(displayed):
    tool = two_events_tool(onset_selection=OnsetSelection.INTERACTIVE)
    tool.run_data()
    box = tool.run_onsets()
    wgt_event = find(box, widgets.Dropdown, "Event:")[0]
    wgt_channel = find(box, widgets.Dropdown, "Channel:")[0]
    wgt_viewing = find(box, widgets.Dropdown, "Selected viewing:")[0]
    btn_previous, btn_next = find(box, description="◀ Previous")[0], find(box, description="Next ▶")[0]
    wgt_status = box.children[1].children[1]
    selected = tool.parameters.selected_onsets

    # the first channel of the first event, with the viewings with an onset
    assert [label for label, _ in wgt_channel.options] == ["HET protons 2.00-8.00 MeV", "HET electrons 0.50-2.00 MeV"]
    assert wgt_viewing.options == (("None", None), ("sun", "sun"))
    assert wgt_viewing.value == "sun"
    assert wgt_status.value == "Event 1: 2 of 2 channels with a selected viewing | channel 1 of 4"
    assert btn_previous.disabled and not btn_next.disabled

    wgt_viewing.value = None
    assert selected.loc[wgt_channel.value, "Viewing"] is None
    assert wgt_status.value.startswith("Event 1: 1 of 2 channels")

    btn_next.click()
    btn_next.click()
    # the first channel of the second event
    assert wgt_event.value == 2
    assert wgt_status.value.endswith("channel 3 of 4")
    btn_next.click()
    assert btn_next.disabled
    btn_previous.click()
    assert wgt_status.value.endswith("channel 3 of 4")

    wgt_event.value = 1
    assert wgt_status.value.endswith("channel 1 of 4")
    wgt_channel.value = wgt_channel.options[1][1]
    assert wgt_status.value.endswith("channel 2 of 4")
    assert selected["Viewing"].tolist() == [None, "sun", "sun", "sun"]


def test_run_vda(displayed, monkeypatch, tmp_path, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(VDA, "define_spacecraft_parameters", lambda self: None)
    monkeypatch.setattr(observers.SolarOrbiter, "distance", lambda self, time: 0.8 * u.AU)
    tool = two_events_tool()
    tool.run_data()
    assert tool.run_onsets() is None
    tool.run_vda()
    assert "Event 1 (2021-10-28 14:00 to 2021-10-28 20:00)" in capsys.readouterr().out
    # the tables of the steps, and one VDA plot per event, saved
    assert [type(o).__name__ for o in displayed] == ["DataFrame", "DataFrame", "DataFrame", "Image", "Image", "DataFrame"]
    assert sorted(f.name for f in tmp_path.glob("*.png")) == [
        "2021-10-28_1400_2021-10-28_2000_het-protons_het-electrons_5min.png",
        "2021-11-09_1500_2021-11-09_2100_het-protons_het-electrons_5min.png",
    ]
