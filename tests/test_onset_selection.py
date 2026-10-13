"""Onset selection"""
import pandas as pd
import pytest

from tests.helpers import make_vda


def vda_with_onset_options(viewings):
    """Two channels: onsets in sun and north for the first, only in north for the second"""
    v = make_vda(viewings=viewings)
    rows = [(1, "het", "protons", "H_Flux", "HET/protons Channel 1", "sun"),
            (1, "het", "protons", "H_Flux", "HET/protons Channel 1", "north"),
            (1, "het", "protons", "H_Flux", "HET/protons Channel 2", "north")]
    v.df_options = pd.DataFrame(
        {"Onset Time": [pd.Timestamp("2021-10-28 16:00")] * len(rows)},
        index=pd.MultiIndex.from_tuples(rows, names=["Event No", "sensor", "particle", "prefix", "channels", "viewing"]),
    )
    return v


@pytest.mark.parametrize("viewings, expected", [
    (["sun", "north"], ["sun", "north"]),
    (["north", "sun"], ["north", "north"]),
])
def test_use_all_follows_viewings_order(viewings, expected):
    v = vda_with_onset_options(viewings)
    v.select_onsets()
    assert v.parameters.selected_onsets["Viewing"].tolist() == expected
