"""Grouping of the energy channels, saved data and the folder of the downloaded data"""
from datetime import datetime

import pandas as pd
import pytest

from spearhead.vda import VDA, VDA_parameters
from spearhead.vda import observers
from tests.helpers import make_vda, rising_epd_data


def test_group_channels_weighted_by_bin_width():
    index = pd.date_range("2021-10-28 14:00", periods=3, freq="5min")
    df = pd.DataFrame({"H_Flux_0": [1.0, 2.0, 3.0], "H_Flux_1": [3.0, 4.0, 5.0], "H_Flux_2": [10.0, 10.0, 10.0]}, index=index)
    grouped = VDA(VDA_parameters())._group_channels_de(
        df, [["H_Flux_0", "H_Flux_1"], ["H_Flux_2"]], [[1.0, 3.0], [2.0]], ["group 1", "group 2"]
    )
    # I = sum(I_n * dE_n) / sum(dE_n)
    assert grouped["group 1"].tolist() == pytest.approx([(1 * 1 + 3 * 3) / 4, (2 * 1 + 4 * 3) / 4, (3 * 1 + 5 * 3) / 4])
    assert grouped["group 2"].tolist() == pytest.approx([10.0, 10.0, 10.0])


def test_data_times_in_order_with_a_gap(monkeypatch):
    # the sun viewing starts 10 minutes later: its times come first, but the data are in time order
    def epd_load(**kwargs):
        data = rising_epd_data(**kwargs)
        if kwargs["viewing"] == "sun":
            data = tuple(df[df.index >= "2021-10-28 14:10"] for df in data[:2]) + (data[2],)
        return data

    monkeypatch.setattr(observers, "epd_load", epd_load)
    v = make_vda(date_ranges=[(datetime(2021, 10, 28, 14), datetime(2021, 10, 28, 15))], resample_frequency="",
                 viewings=["sun", "north"],
                 channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}}})
    v.construct_times_df()
    v.construct_energies_df()
    v.construct_particles_df()
    times = v.df_data.loc[1].index
    assert times.is_monotonic_increasing
    assert times[0] == pd.Timestamp("2021-10-28 14:00") and len(times) == 61
    assert v.df_data.loc[1]["het"]["protons"]["sun"]["H_Flux"].loc[:"2021-10-28 14:09"].isna().all().all()


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


@pytest.mark.parametrize("data_path, expected", [
    ("particle_data", "{cwd}/particle_data"),
    ("", "{cwd}/particle_data"),
    ("data/sep", "{cwd}/data/sep"),
    ("{tmp}/shared", "{tmp}/shared"),
    ("~/sep_data", "{home}/sep_data"),
])
def test_data_path(tmp_path, monkeypatch, data_path, expected):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    v = make_vda(data_path=data_path.format(tmp=tmp_path))
    assert v.DATA_PATH == expected.format(cwd=tmp_path, tmp=tmp_path, home=tmp_path / "home")


def test_data_path_relative_to_the_current_directory_of_the_download(tmp_path, monkeypatch):
    paths = []

    def load(self, sensor, startdate, enddate, viewing, path, particles=None):
        paths.append(path)
        return {}

    monkeypatch.setattr(observers.StereoA, "load", load)
    v = make_vda(observer="sta", data_path="sep")
    monkeypatch.chdir(tmp_path)
    v._load_observer_data("het", "2021-10-28 14:00", "2021-10-28 15:00", "omni", ["protons"])
    assert paths == [str(tmp_path / "sep")]
