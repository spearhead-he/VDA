"""Observers (spacecraft) of the VDA analysis.

Each observer describes its sensors, particles, energy channels and viewings, loads its data in a common
format, and gives its distance from the Sun for the light travel time:
- load(): flux of each particle with the channel numbers as columns, and the energy bins of the channels
- distance(): Sun-observer distance at a given time, after initialize_position()

Sensors without viewings have the single viewing "omni".

OBSERVERS holds the available observers by name, the values of VDA_parameters.observer.
"""
import os
import re
from typing import NamedTuple

import numpy as np
import pandas as pd
from seppy.loader.psp import psp_isois_load
from seppy.loader.soho import soho_load
from seppy.loader.stereo import stereo_load
from seppy.loader.wind import wind3dp_load
from sunpy.coordinates import get_horizons_coord, spice
from sunpy.data import cache
from solo_epd_loader import epd_load


class ParticleData(NamedTuple):
    # Flux with the channel numbers as columns
    flux: pd.DataFrame
    # "Low Energy" and "Bin Width" (MeV) of each channel number
    energies: pd.DataFrame


# Number of channels of the default grouped channels of each particle
DEFAULT_GROUP_SIZE = {"protons": 3, "electrons": 2}


def default_groups(channels, particle) -> list:
    """Consecutive channels in groups of DEFAULT_GROUP_SIZE. A single channel left at the end joins the previous group"""
    size = DEFAULT_GROUP_SIZE[particle]
    groups = [list(channels[i:i + size]) for i in range(0, len(channels), size)]
    if len(groups) > 1 and len(groups[-1]) == 1:
        last = groups.pop()
        groups[-1] += last
    return groups


class Observer:
    # Name used in VDA_parameters.observer
    name: str
    label: str
    # Particles of each sensor
    SENSORS_PARTICLES: dict
    # Channel numbers of each sensor and particle
    CHANNELS: dict
    # All the viewings, and the viewings of each sensor
    VIEWINGS: tuple
    SENSOR_VIEWINGS: dict
    DEFAULT_VIEWINGS: tuple
    # Descriptions of viewings whose names do not say their direction, shown in the parameters form
    VIEWING_DESCRIPTIONS: dict = {}
    # Default grouped channels: {particle: {sensor: [[channel numbers], ...]}}
    DEFAULT_CHANNEL_GROUPS: dict
    # Name of the observer in JPL Horizons, for the default distance()
    HORIZONS_NAME: str | None = None
    # Dates when the energy ranges of the channels of a sensor changed: {sensor: [(date, reason), ...]}.
    # The energy ranges are read from the first event, so all the events must be on the same side of each date
    ENERGY_CHANGES: dict = {}

    def load(self, sensor, startdate, enddate, viewing, path, particles=None) -> dict[str, ParticleData]:
        """Returns the data of the particles of the sensor (by default all of them)"""
        raise NotImplementedError

    def initialize_position(self) -> None:
        """Prepares distance(), e.g. by loading the SPICE kernels"""

    def distance(self, time):
        """Returns the Sun-observer distance (astropy Quantity) at the time, by default from JPL Horizons"""
        return get_horizons_coord(self.HORIZONS_NAME, time).radius


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
    SENSOR_VIEWINGS = {"het": VIEWINGS, "ept": VIEWINGS}
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

    def load(self, sensor, startdate, enddate, viewing, path, particles=None) -> dict[str, ParticleData]:
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
            if particle not in (particles or self.SENSORS_PARTICLES[sensor]):
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


class StereoA(Observer):
    name = "sta"
    label = "STEREO-A"

    # the "protons" of SEPT are ions
    SENSORS_PARTICLES = {
        "het": ("protons", "electrons"),
        "sept": ("protons", "electrons"),
    }

    CHANNELS = {
        "het": {
            "protons": tuple(range(11)),
            "electrons": tuple(range(3)),
        },
        "sept": {
            "protons": tuple(range(2, 32)),
            "electrons": tuple(range(2, 17)),
        },
    }

    VIEWINGS = ("sun", "asun", "north", "south", "omni")
    SENSOR_VIEWINGS = {"het": ("omni",), "sept": ("sun", "asun", "north", "south")}
    DEFAULT_VIEWINGS = ("sun", "omni")
    HORIZONS_NAME = "STEREO-A"

    DEFAULT_CHANNEL_GROUPS = {
        "protons": {"het": default_groups(CHANNELS["het"]["protons"], "protons")},
        "electrons": {"het": default_groups(CHANNELS["het"]["electrons"], "electrons")},
    }

    # Species of the STEREO loader per particle, and flux columns of HET
    SPECIES = {"protons": "p", "electrons": "e"}
    HET_FLUX_COLUMN = {"protons": "Proton_Flux", "electrons": "Electron_Flux"}

    @staticmethod
    def _particle_data(flux, channels_dict_df) -> ParticleData:
        # e.g. "ch_12" (SEPT) or "Proton_Flux_12" (HET) is channel 12
        flux = flux.set_axis([int(c.rsplit("_", 1)[1]) for c in flux.columns], axis="columns")
        return ParticleData(
            flux,
            pd.DataFrame(
                {
                    "Low Energy": channels_dict_df.loc[flux.columns, "lower_E"].to_numpy(),
                    "Bin Width": channels_dict_df.loc[flux.columns, "DE"].to_numpy(),
                },
                index=flux.columns,
            ),
        )

    def load(self, sensor, startdate, enddate, viewing, path, particles=None) -> dict[str, ParticleData]:
        path = os.path.join(path, "stereo")
        os.makedirs(path, exist_ok=True)
        particles = particles or self.SENSORS_PARTICLES[sensor]

        def no_data():
            return ValueError(f"No {self.label} {sensor.upper()} data from {startdate} to {enddate}")

        data = {}
        if sensor == "het":
            df, meta = stereo_load("HET", startdate, enddate, spacecraft="ahead", path=path)
            if len(df) == 0:
                raise no_data()
            for particle in particles:
                columns = [c for c in df.columns if c.startswith(f"{self.HET_FLUX_COLUMN[particle]}_")]
                data[particle] = self._particle_data(df[columns], meta[f"channels_dict_df_{self.SPECIES[particle]}"])
        else:
            # SEPT has one file per species
            for particle in particles:
                species = self.SPECIES[particle]
                df, meta = stereo_load("SEPT", startdate, enddate, spacecraft="ahead",
                                       sept_species=species, sept_viewing=viewing, path=path)
                if len(df) == 0:
                    raise no_data()
                columns = [c for c in df.columns if c.startswith("ch_")]
                data[particle] = self._particle_data(df[columns], meta[f"channels_dict_df_{species}"])
        return data


class ParkerSolarProbe(Observer):
    name = "psp"
    label = "Parker Solar Probe"

    # IMPACT/EPI-Hi HET
    SENSORS_PARTICLES = {"het": ("protons", "electrons")}

    # the proton channels 0-2 and 12-14 have no data
    CHANNELS = {
        "het": {
            "protons": tuple(range(3, 12)),
            "electrons": tuple(range(19)),
        },
    }

    # apertures A and B: the particle flow directions of the data (HET_A_RTN, HET_B_RTN) are about +R and -R,
    # close to the nominal Parker spiral (HET_A_SA, HET_B_SA)
    VIEWINGS = ("A", "B")
    SENSOR_VIEWINGS = {"het": VIEWINGS}
    DEFAULT_VIEWINGS = ("A",)
    VIEWING_DESCRIPTIONS = {
        "A": "sunward, as the sun viewing",
        "B": "anti-sunward, as the asun viewing",
    }
    HORIZONS_NAME = "Parker Solar Probe"

    DEFAULT_CHANNEL_GROUPS = {
        "protons": {"het": default_groups(CHANNELS["het"]["protons"], "protons")},
        "electrons": {"het": default_groups(CHANNELS["het"]["electrons"], "electrons")},
    }

    DATASET = "PSP_ISOIS-EPIHI_L2-HET-RATES60"
    # Energy keys and data columns of each particle. The electrons have count rates, not fluxes
    ENERGY_KEY = {"protons": "H", "electrons": "Electrons"}
    DATA_COLUMN = {"protons": "H_Flux", "electrons": "Electrons_Rate"}

    def load(self, sensor, startdate, enddate, viewing, path, particles=None) -> dict[str, ParticleData]:
        path = os.path.join(path, "psp")
        os.makedirs(path, exist_ok=True)
        df, energies = psp_isois_load(self.DATASET, startdate, enddate, path=path)
        if len(df) == 0:
            raise ValueError(f"No {self.label} EPI-Hi HET data from {startdate} to {enddate}")
        data = {}
        for particle in particles or self.SENSORS_PARTICLES[sensor]:
            channels = list(self.CHANNELS[sensor][particle])
            key = self.ENERGY_KEY[particle]
            delta_minus = np.asarray(energies[f"{key}_ENERGY_DELTAMINUS"])[channels]
            width = delta_minus + np.asarray(energies[f"{key}_ENERGY_DELTAPLUS"])[channels]
            flux = df[[f"{viewing}_{self.DATA_COLUMN[particle]}_{c}" for c in channels]].set_axis(channels, axis="columns")
            if particle == "electrons":
                # count rates per MeV, proportional to the intensity
                flux = flux / width
            data[particle] = ParticleData(
                flux,
                pd.DataFrame(
                    {
                        "Low Energy": np.asarray(energies[f"{key}_ENERGY"])[channels] - delta_minus,
                        "Bin Width": width,
                    },
                    index=channels,
                ),
            )
        return data


class Soho(Observer):
    name = "soho"
    label = "SOHO"

    # ERNE-HED protons and COSTEP-EPHIN electrons
    SENSORS_PARTICLES = {
        "erne": ("protons",),
        "ephin": ("electrons",),
    }

    # EPHIN electron channels: 0 is E150 and 2 is E1300. E300 (1) is deactivated since the failure mode D
    # (4 Oct 2017) and E3000 (3) has no data
    CHANNELS = {
        "erne": {"protons": tuple(range(7))},
        "ephin": {"electrons": (0, 2)},
    }

    VIEWINGS = ("omni",)
    SENSOR_VIEWINGS = {"erne": VIEWINGS, "ephin": VIEWINGS}
    DEFAULT_VIEWINGS = ("omni",)
    HORIZONS_NAME = "SOHO"
    # E1300: 2.64-10.4 MeV before, 0.67-10.4 MeV since
    ENERGY_CHANGES = {"ephin": [(pd.Timestamp("2017-10-04"), "the failure mode D of EPHIN")]}

    DEFAULT_CHANNEL_GROUPS = {
        "protons": {"erne": default_groups(CHANNELS["erne"]["protons"], "protons")},
        "electrons": {"ephin": default_groups(CHANNELS["ephin"]["electrons"], "electrons")},
    }

    EPHIN_ELECTRON_COLUMNS = ("E150", "E300", "E1300", "E3000")

    def load(self, sensor, startdate, enddate, viewing, path, particles=None) -> dict[str, ParticleData]:
        path = os.path.join(path, "soho")
        os.makedirs(path, exist_ok=True)
        if sensor == "erne":
            df, meta = soho_load("SOHO_ERNE-HED_L2-1MIN", startdate, enddate, path=path)
        else:
            df, meta = soho_load("SOHO_COSTEP-EPHIN_L2-1MIN", startdate, enddate, path=path)
        if len(df) == 0:
            raise ValueError(f"No {self.label} {sensor.upper()} data from {startdate} to {enddate}")

        channels = list(self.CHANNELS[sensor][self.SENSORS_PARTICLES[sensor][0]])
        if sensor == "erne":
            flux = df[[f"PH_{c}" for c in channels]].set_axis(channels, axis="columns")
            energies = meta["channels_dict_df_p"].loc[channels]
            low_energy, width = energies["lower_E"].to_numpy(), energies["DE"].to_numpy()
            return {"protons": ParticleData(
                flux, pd.DataFrame({"Low Energy": low_energy, "Bin Width": width}, index=channels))}

        columns = [self.EPHIN_ELECTRON_COLUMNS[c] for c in channels]
        flux = df[columns].set_axis(channels, axis="columns")
        # the energy ranges of the channels depend on the date, e.g. "0.67 - 10.4 MeV"
        ranges = [[float(e) for e in re.findall(r"\d+(?:\.\d+)?", meta["energy_labels"][c])] for c in columns]
        return {"electrons": ParticleData(
            flux,
            pd.DataFrame({"Low Energy": [low for low, _ in ranges], "Bin Width": [high - low for low, high in ranges]},
                         index=channels),
        )}


class Wind(Observer):
    name = "wind"
    label = "Wind"

    # 3DP solid state telescopes: SST Open (protons) and SST Foil (electrons), omnidirectional fluxes
    SENSORS_PARTICLES = {"3dp": ("protons", "electrons")}

    CHANNELS = {
        "3dp": {
            "protons": tuple(range(9)),
            "electrons": tuple(range(7)),
        },
    }

    VIEWINGS = ("omni",)
    SENSOR_VIEWINGS = {"3dp": VIEWINGS}
    DEFAULT_VIEWINGS = ("omni",)
    HORIZONS_NAME = "Wind"

    DEFAULT_CHANNEL_GROUPS = {
        "protons": {"3dp": default_groups(CHANNELS["3dp"]["protons"], "protons")},
        "electrons": {"3dp": default_groups(CHANNELS["3dp"]["electrons"], "electrons")},
    }

    DATASET = {"protons": "WI_SOSP_3DP", "electrons": "WI_SFSP_3DP"}

    def load(self, sensor, startdate, enddate, viewing, path, particles=None) -> dict[str, ParticleData]:
        path = os.path.join(path, "wind")
        os.makedirs(path, exist_ok=True)
        data = {}
        for particle in particles or self.SENSORS_PARTICLES[sensor]:
            df, meta = wind3dp_load(self.DATASET[particle], startdate, enddate, resample=None, path=path)
            if len(df) == 0:
                raise ValueError(f"No {self.label} 3DP {particle} data from {startdate} to {enddate}")
            channels = list(self.CHANNELS[sensor][particle])
            # fluxes per eV, as intensities per MeV
            flux = df[[f"FLUX_{c}" for c in channels]].set_axis(channels, axis="columns") * 1e6
            # mean energies of the loaded data, with a width of 60% of the mean energy (as in seppy)
            energies = meta["channels_dict_df"].loc[[f"ENERGY_{c}" for c in channels]]
            data[particle] = ParticleData(
                flux,
                pd.DataFrame(
                    {"Low Energy": energies["lower_E"].to_numpy(), "Bin Width": energies["DE"].to_numpy()},
                    index=channels,
                ),
            )
        return data


OBSERVERS = {observer.name: observer for observer in (SolarOrbiter(), StereoA(), ParkerSolarProbe(), Soho(), Wind())}
DEFAULT_OBSERVER = "solo"
