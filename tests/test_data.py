"""Grouping of the energy channels and saved data"""
import pandas as pd
import pytest

from spearhead.vda import VDA, VDA_parameters
from tests.helpers import make_vda


def test_group_channels_weighted_by_bin_width():
    index = pd.date_range("2021-10-28 14:00", periods=3, freq="5min")
    df = pd.DataFrame({"H_Flux_0": [1.0, 2.0, 3.0], "H_Flux_1": [3.0, 4.0, 5.0], "H_Flux_2": [10.0, 10.0, 10.0]}, index=index)
    grouped = VDA(VDA_parameters())._group_channels_de(
        df, [["H_Flux_0", "H_Flux_1"], ["H_Flux_2"]], [[1.0, 3.0], [2.0]], ["group 1", "group 2"]
    )
    # I = sum(I_n * dE_n) / sum(dE_n)
    assert grouped["group 1"].tolist() == pytest.approx([(1 * 1 + 3 * 3) / 4, (2 * 1 + 4 * 3) / 4, (3 * 1 + 5 * 3) / 4])
    assert grouped["group 2"].tolist() == pytest.approx([10.0, 10.0, 10.0])


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
