"""Parameters"""
import dataclasses
from datetime import datetime

import pytest

from spearhead.vda import OnsetSelection, VDA_parameters
from spearhead.vda.conf import DEFAULT_DATA_PATH, check_resample_frequency
from tests.helpers import make_vda


@pytest.mark.parametrize("name, version", [
    ("input_type", "v0.3.0"), ("date_range_filepath", "v0.3.0"), ("reference_times_filepath", "v0.3.0"),
    ("date_start", "v0.3.0"), ("date_end", "v0.3.0"),
    ("viewings_tt", "v0.4.0"), ("load_data", "v0.4.0"), ("save_data", "v0.4.0"),
    ("default_channel_groups", "v0.5.0"),
])
def test_removed_parameters(name, version):
    with pytest.raises(AttributeError, match=f"removed in {version}"):
        setattr(VDA_parameters(), name, None)


def test_unknown_parameter():
    with pytest.raises(AttributeError, match="no parameter 'input_filepth'"):
        VDA_parameters().input_filepth = "events.csv"


def test_onset_selection_values():
    p = VDA_parameters()
    p.onset_selection = 1
    assert p.onset_selection is OnsetSelection.INTERACTIVE
    with pytest.raises(ValueError, match="Unknown onset selection 2"):
        p.onset_selection = 2


def test_parameters_instances_are_independent():
    p = VDA_parameters()
    q = dataclasses.replace(p, input_filepath="events.csv")
    assert p.input_filepath == "" and q.input_filepath == "events.csv"
    # new parameter sets do not share their lists and dictionaries
    other = VDA_parameters()
    other.date_ranges.append((datetime(2021, 1, 1), datetime(2021, 1, 2)))
    other.channel_groups["protons"]["HET/protons Channel 1"]["channels"].append(4)
    assert len(p.date_ranges) == 1 and p.channel_groups["protons"]["HET/protons Channel 1"]["channels"] == [1, 2, 3]
    assert VDA_parameters() == VDA_parameters()


def test_onset_method_parameters_follow_the_method():
    p = VDA_parameters()
    assert p.onset_method_parameters == {"s": 3, "n": 3}
    p.onset_method_parameters["s"] = 4
    p.onset_method = "sigma"
    assert p.onset_method_parameters == {"s": 4, "n": 3}
    p.onset_method = "poisson_cusum"
    assert p.onset_method_parameters == {"cusum_minutes": 30, "sigma_multiplier": 2}
    assert VDA_parameters(onset_method="poisson_cusum").onset_method_parameters == {"cusum_minutes": 30, "sigma_multiplier": 2}
    with pytest.raises(ValueError, match="Unknown onset method 'cusum'"):
        p.onset_method = "cusum"


def test_onset_method_parameters_of_another_method():
    v = make_vda()
    v.parameters.onset_method_parameters = {"cusum_minutes": 30}
    with pytest.raises(ValueError, match="not parameters of the sigma onset method"):
        v.calculate_onsets()


def test_default_channel_groups():
    p = VDA_parameters()
    assert list(p.channel_groups) == ["protons", "electrons"]
    assert len(p.channel_groups["protons"]) == 8
    assert p.channel_groups["electrons"] == {
        "HET/electrons Channel 1": {"sensor": "het", "channels": [0, 1]},
        "HET/electrons Channel 2": {"sensor": "het", "channels": [2, 3]},
    }
    assert p.sensors_particles == {"het": ["protons", "electrons"]}
    assert p.viewings == ["sun"]


def test_default_data_path():
    assert VDA_parameters().data_path == DEFAULT_DATA_PATH == "particle_data"


@pytest.mark.parametrize("frequency", ["", "5min", "30s", "1h"])
def test_valid_resample_frequency(frequency):
    check_resample_frequency(frequency)


@pytest.mark.parametrize("frequency", ["5 minutes", "5", "five"])
def test_invalid_resample_frequency(frequency):
    with pytest.raises(ValueError, match=f"Invalid resample frequency '{frequency}': use a pandas offset alias "
                                         r"\(e.g. 5min, 30s or 1h\), or leave it blank for no resampling"):
        check_resample_frequency(frequency)


def test_intermediate_tables_not_displayed_by_default():
    assert VDA_parameters().view_dfs is False
