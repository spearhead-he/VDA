"""Observers (spacecraft) of the VDA analysis.

Each observer describes its sensors, particles, energy channels and viewings, loads its data in a common
format, and gives its distance from the Sun for the light travel time:
- load(): flux of each particle with the channel numbers as columns, and the energy bins of the channels
- distance(): Sun-observer distance at a given time, after initialize_position()

OBSERVERS holds the available observers by name, the values of VDA_parameters.observer.
"""
from typing import NamedTuple

import pandas as pd
from sunpy.coordinates import spice
from sunpy.data import cache
from solo_epd_loader import epd_load


class ParticleData(NamedTuple):
    # Flux with the channel numbers as columns
    flux: pd.DataFrame
    # "Low Energy" and "Bin Width" (MeV) of each channel number
    energies: pd.DataFrame


class Observer:
    # Name used in VDA_parameters.observer
    name: str
    label: str
    # Particles of each sensor
    SENSORS_PARTICLES: dict
    # Channel numbers of each sensor and particle
    CHANNELS: dict
    VIEWINGS: tuple
    DEFAULT_VIEWINGS: tuple
    # Default grouped channels: {particle: {sensor: [[channel numbers], ...]}}
    DEFAULT_CHANNEL_GROUPS: dict

    def load(self, sensor, startdate, enddate, viewing, path) -> dict[str, ParticleData]:
        """Returns the data of each particle of the sensor"""
        raise NotImplementedError

    def initialize_position(self) -> None:
        """Prepares distance(), e.g. by loading the SPICE kernels"""

    def distance(self, time):
        """Returns the Sun-observer distance (astropy Quantity) at the time"""
        raise NotImplementedError


class SolarOrbiter(Observer):
    name = "solo"
    label = "Solar Orbiter"

    SENSORS_PARTICLES = {
        "het": ("protons", "electrons"),
        "ept": ("protons", "electrons"),
    }

    CHANNELS = {
        "het": {
            "protons": tuple(range(36)),
            "electrons": tuple(range(4)),
        },
        "ept": {
            "protons": tuple(range(64)),
            "electrons": tuple(range(34)),
        },
    }

    VIEWINGS = ("sun", "asun", "north", "south", "omni")
    DEFAULT_VIEWINGS = ("sun",)

    DEFAULT_CHANNEL_GROUPS = {
        "protons": {
            "het": [
                [1, 2, 3],
                [10, 11, 12],
                [13, 14, 15],
                [16, 17, 18],
                [19, 20, 21],
                [22, 23, 24],
                [25, 26, 27],
                [28, 29, 30, 31],
            ],
        },
        "electrons": {
            "het": [
                [0, 1],
                [2, 3],
            ],
        },
    }

    # Flux columns of the EPD data per sensor and particle
    FLUX_COLUMN = {
        "het": {"protons": "H_Flux", "electrons": "Electron_Flux"},
        "ept": {"protons": "Ion_Flux", "electrons": "Electron_Flux"},
    }

    # Energy bins keys of the EPD data per sensor and particle
    ENERGY_BINS_COLUMN = {
        "het": {"protons": "H_Bins", "electrons": "Electron_Bins"},
        "ept": {"protons": "Ion_Bins", "electrons": "Electron_Bins"},
    }

    SPICE_KERNELS_URL = "https://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels"
    SPICE_KERNELS = (
        "ck/solo_ANC_soc-sc-fof-ck_20180930-21000101_V03.bc",
        "ck/solo_ANC_soc-stix-ck_20180930-21000101_V03.bc",
        "ck/solo_ANC_soc-flown-att_20221011T142135-20221012T141817_V01.bc",
        "fk/solo_ANC_soc-sc-fk_V09.tf",
        "fk/solo_ANC_soc-sci-fk_V08.tf",
        "ik/solo_ANC_soc-stix-ik_V02.ti",
        "lsk/naif0012.tls",
        "pck/pck00010.tpc",
        "sclk/solo_ANC_soc-sclk_20231015_V01.tsc",
        "spk/de421.bsp",
        "spk/solo_ANC_soc-orbit-stp_20200210-20301120_280_V1_00288_V01.bsp",
    )

    def load(self, sensor, startdate, enddate, viewing, path) -> dict[str, ParticleData]:
        df_protons, df_electrons, energies = epd_load(
            sensor=sensor,
            level="l2",
            startdate=startdate,
            enddate=enddate,
            viewing=viewing,
            path=path,
            autodownload=True,
        )
        data = {}
        for particle, df_particle in (("protons", df_protons), ("electrons", df_electrons)):
            if particle not in self.SENSORS_PARTICLES[sensor]:
                continue
            flux = df_particle[self.FLUX_COLUMN[sensor][particle]]
            # e.g. "H_Flux_12" is channel 12
            flux = flux.set_axis([int(c.rsplit("_", 1)[1]) for c in flux.columns], axis="columns")
            energy_bins = self.ENERGY_BINS_COLUMN[sensor][particle]
            data[particle] = ParticleData(
                flux,
                pd.DataFrame(
                    {
                        "Low Energy": energies[f"{energy_bins}_Low_Energy"],
                        "Bin Width": energies[f"{energy_bins}_Width"],
                    },
                    index=flux.columns,
                ),
            )
        return data

    def initialize_position(self) -> None:
        spice.initialize([cache.download(f"{self.SPICE_KERNELS_URL}/{kernel}") for kernel in self.SPICE_KERNELS])

    def distance(self, time):
        return spice.get_body("Solar Orbiter", time, spice_frame="SOLO_HEEQ").distance


OBSERVERS = {observer.name: observer for observer in (SolarOrbiter(),)}
DEFAULT_OBSERVER = "solo"
