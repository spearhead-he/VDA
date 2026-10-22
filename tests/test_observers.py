"""Observers (spacecraft): data loading, viewings and default grouped channels"""
import io
import logging
from datetime import datetime

import numpy as np
import pandas as pd
import pytest

from spearhead.vda import VDA_parameters
from spearhead.vda import observers
from tests.helpers import epd_like_data, make_vda


def test_solar_orbiter_load(monkeypatch, tmp_path):
    index = pd.date_range("2021-10-28 14:00", periods=4, freq="1min")
    paths = []

    def epd_load(**kwargs):
        paths.append(kwargs["path"])
        return epd_like_data(index, kwargs["sensor"])

    monkeypatch.setattr(observers, "epd_load", epd_load)
    data = observers.SolarOrbiter().load("ept", index[0], index[-1], "sun", str(tmp_path))
    # the files of each observer are in its subfolder
    assert paths == [str(tmp_path / "solo")]
    assert (tmp_path / "solo").is_dir()
    assert list(data) == ["protons", "electrons"]
    assert data["protons"].flux.columns.tolist() == [0, 1, 2]
    assert data["protons"].flux[2].tolist() == epd_like_data(index, "ept")[0][("Ion_Flux", "Ion_Flux_2")].tolist()
    assert data["electrons"].energies.loc[1].tolist() == [1.0, 1.0]


def test_data_from_the_observer(monkeypatch):
    # the observer data become the columns of the analysis, e.g. H_Flux_2 for proton channel 2
    index = pd.date_range("2021-10-28 14:00", periods=20, freq="1min")
    monkeypatch.setattr(observers, "epd_load", lambda **kwargs: epd_like_data(index, kwargs["sensor"]))
    v = make_vda(date_ranges=[(index[0], index[-1])], resample_frequency="5min", viewings=["sun", "north"],
                 channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [1, 2]}}})
    v.construct_times_df()
    v.construct_energies_df()
    assert v.df_energies.loc[("het", "H_Flux_2")].tolist() == [4.0, 4.0, 8.0]
    assert v.df_energies.loc[("ept", "Electron_Flux_1"), "High Energy"] == 2.0
    v.construct_particles_df()
    assert v.df_data.attrs["observer"] == "solo"
    assert v.df_data.columns.tolist() == [("het", "protons", viewing, "H_Flux", f"H_Flux_{c}")
                                          for viewing in ("sun", "north") for c in range(3)]
    assert len(v.df_data) == 4


def stereo_like_data(instrument, startdate, enddate, spacecraft, sept_species=None, sept_viewing=None, path=None):
    """Data as returned by stereo_load: HET with 2 proton and 2 electron channels, SEPT with channels 2 to 4"""
    index = pd.date_range("2021-10-28 14:00", periods=20, freq="1min")

    def channels_dict_df(channels):
        lower = np.array([1.0, 2.0, 4.0])[:len(channels)]
        return pd.DataFrame({"lower_E": lower, "upper_E": 2 * lower, "DE": lower}, index=channels)

    if instrument == "HET":
        columns = ["Electron_Flux_0", "Electron_Flux_1", "Electron_Sigma_0", "Electron_Sigma_1",
                   "Proton_Flux_0", "Proton_Flux_1", "Proton_Sigma_0", "Proton_Sigma_1"]
        meta = {"channels_dict_df_e": channels_dict_df([0, 1]), "channels_dict_df_p": channels_dict_df([0, 1])}
    else:
        assert sept_viewing in ("sun", "asun", "north", "south")
        columns = ["ch_2", "ch_3", "ch_4", "err_ch_2", "err_ch_3", "err_ch_4"]
        meta = {f"channels_dict_df_{sept_species}": channels_dict_df([2, 3, 4])}
    flux = pd.DataFrame(np.ones((len(index), len(columns))), index=index, columns=columns)
    return flux, meta


def test_stereo_a_load(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "stereo_load", stereo_like_data)
    sta = observers.StereoA()
    het = sta.load("het", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path))
    assert list(het) == ["protons", "electrons"]
    assert het["protons"].flux.columns.tolist() == [0, 1]
    assert het["electrons"].energies.loc[1].tolist() == [2.0, 2.0]
    sept = sta.load("sept", "2021-10-28 14:00", "2021-10-28 15:00", "north", str(tmp_path), particles=["electrons"])
    assert list(sept) == ["electrons"]
    assert sept["electrons"].flux.columns.tolist() == [2, 3, 4]
    assert sept["electrons"].energies.loc[4].tolist() == [4.0, 4.0]


def test_stereo_a_without_data(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "stereo_load", lambda *args, **kwargs: ([], []))
    with pytest.raises(ValueError, match="No STEREO-A HET data"):
        observers.StereoA().load("het", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path))


def psp_like_data(dataset, startdate, enddate, path=None):
    """Data as returned by psp_isois_load for EPI-Hi HET: 15 proton and 19 electron channels per aperture"""
    index = pd.date_range("2021-10-28 14:00", periods=20, freq="1min")
    columns = {f"{aperture}_{column}_{c}": 10.0 if aperture == "A" else 20.0
               for aperture in ("A", "B")
               for column, n_channels in (("H_Flux", 15), ("Electrons_Rate", 19))
               for c in range(n_channels)}
    energies = {}
    for key, n_channels in (("H", 15), ("Electrons", 19)):
        energies[f"{key}_ENERGY"] = np.arange(1.0, n_channels + 1)
        energies[f"{key}_ENERGY_DELTAMINUS"] = np.full(n_channels, 0.25)
        energies[f"{key}_ENERGY_DELTAPLUS"] = np.full(n_channels, 0.25)
    return pd.DataFrame(columns, index=index), energies


def test_parker_solar_probe_load(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "psp_isois_load", psp_like_data)
    data = observers.ParkerSolarProbe().load("het", "2021-10-28 14:00", "2021-10-28 15:00", "B", str(tmp_path))
    # the proton channels without data are left out
    assert data["protons"].flux.columns.tolist() == list(range(3, 12))
    assert data["protons"].energies.loc[3].tolist() == [3.75, 0.5]
    assert (data["protons"].flux == 20.0).all().all()
    # the electron count rates are divided by the bin width
    assert data["electrons"].flux.columns.tolist() == list(range(19))
    assert (data["electrons"].flux == 40.0).all().all()


def test_parker_solar_probe_default_groups():
    p = VDA_parameters(observer="psp")
    assert p.viewings == ["A"]
    assert [g["channels"] for g in p.channel_groups["protons"].values()] == [[3, 4, 5], [6, 7, 8], [9, 10, 11]]
    assert [g["channels"] for g in p.channel_groups["electrons"].values()][-1] == [16, 17, 18]


def soho_like_data(dataset, startdate, enddate, path=None):
    """Data as returned by soho_load: ERNE-HED with 7 proton channels, EPHIN after the failure mode D"""
    index = pd.date_range("2021-10-28 14:00", periods=20, freq="1min")
    if dataset == "SOHO_ERNE-HED_L2-1MIN":
        columns = [f"{p}_{c}" for p in ("AH", "PH") for c in range(7)]
        lower = np.array([13.0, 16, 20, 25, 32, 40, 50])
        meta = {"channels_dict_df_p": pd.DataFrame({"lower_E": lower, "DE": np.diff(np.append(lower, 64.0))})}
    else:
        columns = ["E150", "E300", "E1300", "E3000", "P4"]
        meta = {"energy_labels": {"E150": "0.25 - 0.7 MeV", "E300": "deactivated bc. of failure mode D",
                                  "E1300": "0.67 - 10.4 MeV", "E3000": "4.80 - 10.4 MeV", "P4": "4.3 - 7.8 MeV"}}
    return pd.DataFrame({c: float(i) for i, c in enumerate(columns)}, index=index), meta


def test_soho_load(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "soho_load", soho_like_data)
    soho = observers.Soho()
    erne = soho.load("erne", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path))
    assert list(erne) == ["protons"]
    assert erne["protons"].flux.columns.tolist() == list(range(7))
    assert (erne["protons"].flux[0] == 7.0).all()
    assert erne["protons"].energies.loc[6].tolist() == [50.0, 14.0]
    ephin = soho.load("ephin", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path))
    assert list(ephin) == ["electrons"]
    # E150 and E1300, with the energy ranges of their labels
    assert ephin["electrons"].flux.columns.tolist() == [0, 2]
    assert (ephin["electrons"].flux[2] == 2.0).all()
    assert ephin["electrons"].energies["Low Energy"].tolist() == [0.25, 0.67]
    assert ephin["electrons"].energies["Bin Width"].tolist() == pytest.approx([0.45, 9.73])


def test_soho_default_groups():
    p = VDA_parameters(observer="soho")
    assert p.viewings == ["omni"]
    assert p.sensors_particles == {"erne": ["protons"], "ephin": ["electrons"]}
    assert [g["channels"] for g in p.channel_groups["protons"].values()] == [[0, 1, 2], [3, 4, 5, 6]]
    assert [g["channels"] for g in p.channel_groups["electrons"].values()] == [[0, 2]]


def wind_like_data(dataset, startdate, enddate, resample=None, path=None):
    """Data as returned by wind3dp_load for the omnidirectional fluxes: 9 proton or 7 electron channels"""
    n_channels = {"WI_SOSP_3DP": 9, "WI_SFSP_3DP": 7}[dataset]
    index = pd.date_range("2021-10-28 14:00", periods=20, freq="12s")
    columns = {"TIME": 0.0, **{f"FLUX_{c}": 1e-6 * (c + 1) for c in range(n_channels)},
               **{f"ENERGY_{c}": 1e5 * (c + 1) for c in range(n_channels)}}
    mean_e = np.arange(1, n_channels + 1) * 0.1
    energies = pd.DataFrame({"lower_E": 0.7 * mean_e, "DE": 0.6 * mean_e},
                            index=[f"ENERGY_{c}" for c in range(n_channels)])
    return pd.DataFrame(columns, index=index), {"channels_dict_df": energies}


def test_wind_load(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "wind3dp_load", wind_like_data)
    data = observers.Wind().load("3dp", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path))
    assert data["protons"].flux.columns.tolist() == list(range(9))
    assert data["electrons"].flux.columns.tolist() == list(range(7))
    # fluxes per eV become intensities per MeV
    assert data["electrons"].flux[2].tolist() == pytest.approx([3.0] * 20)
    assert data["protons"].energies.loc[1].tolist() == pytest.approx([0.14, 0.12])
    only_electrons = observers.Wind().load("3dp", "2021-10-28 14:00", "2021-10-28 15:00", "omni", str(tmp_path),
                                           particles=["electrons"])
    assert list(only_electrons) == ["electrons"]


def test_wind_default_groups():
    p = VDA_parameters(observer="wind")
    assert p.viewings == ["omni"]
    assert [g["channels"] for g in p.channel_groups["protons"].values()] == [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    assert [g["channels"] for g in p.channel_groups["electrons"].values()] == [[0, 1], [2, 3], [4, 5, 6]]


def bepi_like_data(startdate, enddate=None, path=None):
    """Data as returned by bepi_sixsp_l3_loader: sides 0-3 with P1-P9 and E1-E7, with UTC times"""
    # the loader is called from midnight
    assert startdate == pd.Timestamp("2023-03-13")
    index = pd.date_range("2023-03-13 00:01", periods=20, freq="2min", tz="UTC")
    columns = {f"Side{side}_{p}{c}": float(10 * side + c)
               for side in range(4) for p, n in (("P", 9), ("E", 7)) for c in range(1, n + 1)}
    meta = {}
    for side in range(4):
        for key, p, n in (("Proton", "P", 9), ("Electron", "E", 7)):
            meta[f"Side{side}_{key}_Bins_Low_Energy"] = {f"{p}{c}": float(c) for c in range(1, n + 1)}
            meta[f"Side{side}_{key}_Bins_High_Energy"] = {f"{p}{c}": c + 0.5 + side for c in range(1, n + 1)}
    return pd.DataFrame(columns, index=index), meta


def test_bepicolombo_load(monkeypatch, tmp_path):
    monkeypatch.setattr(observers, "bepi_sixsp_l3_loader", bepi_like_data)
    data = observers.BepiColombo().load("sixs", "2023-03-13 06:00", "2023-03-13 12:00", "side2", str(tmp_path))
    assert data["protons"].flux.columns.tolist() == list(range(1, 10))
    assert data["electrons"].flux.columns.tolist() == list(range(1, 8))
    # the data and energies of the side, with times without time zone
    assert (data["electrons"].flux[3] == 23.0).all()
    assert data["protons"].energies.loc[4].tolist() == [4.0, 2.5]
    assert data["protons"].flux.index.tz is None


def test_bepicolombo_defaults():
    p = VDA_parameters(observer="bepi")
    assert p.viewings == ["side0"]
    assert p.AVAILABLE_VIEWINGS == ("side0", "side1", "side2")
    assert [g["channels"] for g in p.channel_groups["protons"].values()] == [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    assert [g["channels"] for g in p.channel_groups["electrons"].values()] == [[1, 2], [3, 4], [5, 6, 7]]


@pytest.mark.parametrize("date_ranges, channel_groups, error", [
    # before and after the failure mode D of EPHIN
    ([(datetime(2012, 5, 17), datetime(2012, 5, 18)), (datetime(2021, 10, 28), datetime(2021, 10, 29))], None, True),
    # an event across it
    ([(datetime(2017, 10, 3), datetime(2017, 10, 5))], None, True),
    ([(datetime(2021, 10, 28), datetime(2021, 10, 29)), (datetime(2024, 5, 11), datetime(2024, 5, 12))], None, False),
    # without EPHIN channels
    ([(datetime(2012, 5, 17), datetime(2012, 5, 18)), (datetime(2021, 10, 28), datetime(2021, 10, 29))],
     {"protons": {"ERNE/protons Channel 1": {"sensor": "erne", "channels": [0, 1, 2]}}}, False),
])
def test_energy_changes(date_ranges, channel_groups, error):
    v = make_vda(observer="soho", date_ranges=date_ranges)
    if channel_groups:
        v.parameters.channel_groups = channel_groups
    v.construct_times_df()
    if error:
        with pytest.raises(ValueError, match="SOHO EPHIN channels changed on 2017-10-04"):
            v.check_energy_changes()
    else:
        v.check_energy_changes()


def test_viewings_of_each_sensor(monkeypatch):
    # STEREO-A HET has only the omni viewing, SEPT has no omni viewing
    monkeypatch.setattr(observers, "stereo_load", stereo_like_data)
    v = make_vda(observer="sta", viewings=["north", "omni", "sun"],
                 channel_groups={"protons": {"HET/protons Channel 1": {"sensor": "het", "channels": [0, 1]}},
                                 "electrons": {"SEPT/electrons Channel 1": {"sensor": "sept", "channels": [2, 3]}}})
    assert v.sensor_viewings("het") == ["omni"]
    assert v.sensor_viewings("sept") == ["north", "sun"]
    v.construct_times_df()
    v.construct_energies_df()
    v.construct_particles_df()
    assert sorted({c[:3] for c in v.df_data.columns}) == [
        ("het", "protons", "omni"), ("sept", "electrons", "north"), ("sept", "electrons", "sun")]
    v.group_energy_channels()
    v.construct_energy_channels_characteristics()
    assert v.df_channels_chars.index.get_level_values("sensor").tolist() == ["het", "sept"]


def test_sensor_without_selected_viewings():
    v = make_vda(observer="sta", viewings=["sun"])
    with pytest.raises(ValueError, match="None of the selected viewings \\(sun\\) is a viewing of STEREO-A HET"):
        v.check_viewings()


def test_observer_change_resets_viewings_and_channel_groups():
    p = VDA_parameters(viewings=["north"])
    p.observer = "solo"
    assert p.viewings == ["north"]
    p.observer = "sta"
    assert p.viewings == ["sun", "omni"]
    assert p.sensors_particles == {"het": ["protons", "electrons"]}
    assert [g["channels"] for g in p.channel_groups["protons"].values()] == [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10]]
    assert [g["channels"] for g in p.channel_groups["electrons"].values()] == [[0, 1, 2]]
    assert p.AVAILABLE_CHANNELS["sept"]["electrons"] == tuple(range(2, 17))


@pytest.mark.parametrize("n_channels, particle, expected", [
    (6, "protons", [[0, 1, 2], [3, 4, 5]]),
    (8, "protons", [[0, 1, 2], [3, 4, 5], [6, 7]]),
    (7, "protons", [[0, 1, 2], [3, 4, 5, 6]]),
    (4, "electrons", [[0, 1], [2, 3]]),
    (5, "electrons", [[0, 1], [2, 3, 4]]),
    (1, "electrons", [[0]]),
])
def test_default_groups(n_channels, particle, expected):
    # a single channel left at the end joins the previous group
    assert observers.default_groups(tuple(range(n_channels)), particle) == expected


def test_unknown_observer():
    with pytest.raises(ValueError, match="Unknown observer 'other'"):
        VDA_parameters(observer="other")
    p = VDA_parameters()
    with pytest.raises(ValueError, match="Unknown observer 'other'"):
        p.observer = "other"


def test_download_once(monkeypatch, tmp_path, capsys):
    # the files are downloaded once, without progress bars, into the cache
    urls = []

    class Response:
        headers = {}

        def __init__(self, url):
            urls.append(url)

        def __enter__(self):
            return io.BytesIO(b"kernel")

        def __exit__(self, *args):
            return False

    storage = observers.SqliteStorage(str(tmp_path / "cache.db"))
    cache = observers.Cache(observers._QuietDownloader(), storage, tmp_path, expiry=None)
    monkeypatch.setattr(observers, "_cache_storage", storage)
    monkeypatch.setattr(observers, "_cache", cache)
    monkeypatch.setattr(observers, "urlopen", lambda url, timeout: Response(url))
    # the file name of the sunpy cache is read from the response of urllib
    monkeypatch.setattr("sunpy.data.data_manager.cache.urlopen", lambda url: Response(url))
    paths = observers.download_once(["https://example.org/a.tls", "https://example.org/b.bsp"], "test kernels")
    assert capsys.readouterr().out == "Downloading the test kernels (only once)\n"
    assert [p.read_bytes() for p in paths] == [b"kernel", b"kernel"]
    assert observers.download_once(["https://example.org/a.tls", "https://example.org/b.bsp"], "test kernels") == paths
    assert capsys.readouterr().out == ""


def test_spiceypy_logs_only_warnings():
    # set by the package, whatever the level of the root logger (set to INFO by drms)
    assert logging.getLogger("spiceypy").level == logging.WARNING


def test_quiet_downloader_error(monkeypatch, tmp_path):
    def no_connection(url, timeout):
        raise OSError("Connection refused")

    monkeypatch.setattr(observers, "urlopen", no_connection)
    with pytest.raises(observers.DownloaderError, match="Connection refused"):
        observers._QuietDownloader().download("https://example.org/a.tls", str(tmp_path / "a.tls"))
