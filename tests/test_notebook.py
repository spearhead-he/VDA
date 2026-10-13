"""Widgets of the notebook interface"""
import pandas as pd

from spearhead.vda import VDA, VDA_parameters
from spearhead.vda import observers
from spearhead.vda.notebook import VDA_notebook
from tests.helpers import epd_like_data


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
