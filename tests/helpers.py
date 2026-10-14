"""Helpers of the tests: parameters, files and data without downloads, and figures for image comparison"""
import numpy as np
import pandas as pd

from spearhead.vda import VDA, VDA_parameters
from spearhead.vda.observers import ParticleData


def make_vda(**parameters) -> VDA:
    p = VDA_parameters()
    for name, value in parameters.items():
        setattr(p, name, value)
    return VDA(p)


def write_csv(tmp_path, name, text) -> str:
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def epd_like_data(index, sensor="het"):
    """Data as returned by epd_load: protons with 3 channels and electrons with 2"""
    proton_flux, proton_bins = {"het": ("H_Flux", "H_Bins"), "ept": ("Ion_Flux", "Ion_Bins")}[sensor]

    def particle(flux_column, n_channels):
        columns = pd.MultiIndex.from_tuples(
            [(flux_column, f"{flux_column}_{c}") for c in range(n_channels)]
            + [("Uncertainty", f"Uncertainty_{c}") for c in range(n_channels)]
        )
        return pd.DataFrame(np.arange(len(index) * 2 * n_channels, dtype=float).reshape(len(index), -1),
                            index=index, columns=columns)

    energies = {
        f"{proton_bins}_Low_Energy": np.array([1.0, 2.0, 4.0]),
        f"{proton_bins}_Width": np.array([1.0, 2.0, 4.0]),
        "Electron_Bins_Low_Energy": np.array([0.5, 1.0]),
        "Electron_Bins_Width": np.array([0.5, 1.0]),
    }
    return particle(proton_flux, 3), particle("Electron_Flux", 2), energies


def rising_epd_data(**kwargs):
    """Data as returned by epd_load, from the start to the end date at 1 minute, with a rise two hours after
    the start in the sun viewing only"""
    index = pd.date_range(kwargs["startdate"], kwargs["enddate"], freq="1min")
    df_protons, df_electrons, energies = epd_like_data(index, kwargs["sensor"])
    noise = np.where(np.arange(len(index)) % 2 == 0, 1.0, 1.2)
    rise = np.where(index >= index[0] + pd.Timedelta(hours=2), 100.0, 0.0) if kwargs["viewing"] == "sun" else 0.0
    for df in (df_protons, df_electrons):
        for column in df.columns:
            df[column] = noise + rise
    return df_protons, df_electrons, energies


def strip_figure_text(fig):
    """
    Remove all text elements that vary across matplotlib versions.
    Use before returning a figure in mpl_image_compare tests.
    """
    fig.suptitle('')
    for text in fig.texts:
        text.set_text('')
    for ax in fig.axes:
        ax.set_title('')
        ax.set_xlabel('')      # not removed by remove_text=True
        ax.set_ylabel('')      # not removed by remove_text=True
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        leg = ax.get_legend()
        if leg:
            leg.remove()
    return fig


def rising_observer_load(self, sensor, startdate, enddate, viewing, path, particles=None):
    """load() of an observer with every channel of the sensor, from the start to the end date at 1 minute, with a rise
    two hours after the start in every viewing. The proton channel c is from 1 + c MeV, the electron channel c from
    0.1 * (1 + c) MeV"""
    index = pd.date_range(startdate, enddate, freq="1min")
    noise = np.where(np.arange(len(index)) % 2 == 0, 1.0, 1.2)
    rise = np.where(index >= index[0] + pd.Timedelta(hours=2), 100.0, 0.0)
    scale = {"protons": 1.0, "electrons": 0.1}
    data = {}
    for particle in particles or self.SENSORS_PARTICLES[sensor]:
        channels = list(self.CHANNELS[sensor][particle])
        flux = pd.DataFrame({c: noise + rise for c in channels}, index=index)
        energies = pd.DataFrame({"Low Energy": [scale[particle] * (1 + c) for c in channels],
                                 "Bin Width": [scale[particle]] * len(channels)}, index=channels)
        data[particle] = ParticleData(flux, energies)
    return data
