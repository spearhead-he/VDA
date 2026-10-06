import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const

from math import ceil, sqrt
from os import getcwd
from datetime import datetime, timedelta

from matplotlib import pyplot as plt
from matplotlib import dates as mdates
from sunpy.coordinates import spice
from sunpy.data import cache
from solo_epd_loader import epd_load
from pyonset import Onset, BootstrapWindow


class VDA:

    def __init__(self, parameters):
        self.parameters = parameters
        self.results = pd.DataFrame({
            "Release Time": [],
            "Release Time Error": [],
            "Extra Time": [],
            "APL": [],
            "APL Error": []
        })

    ############### Reference Times DF ###############
    EVENT_INDEX_NAME = "Event No"
    REF_TIME_COLNAME = "Reference Time"
    START_TIME_COLNAME = "Start Time"
    BG_START_TIME_COLNAME = "BG Start"
    BG_END_TIME_COLNAME = "BG End"
    END_TIME_COLNAME = "End Time"

    # Accepted columns of the input file, after the event number column
    REFERENCE_TIMES_FILE_LAYOUT = (REF_TIME_COLNAME,)
    DATE_RANGE_FILE_LAYOUTS = (
        # data range only, the default background window is used
        (START_TIME_COLNAME, END_TIME_COLNAME),
        # "Start Time" is also used as the background start
        (START_TIME_COLNAME, BG_END_TIME_COLNAME, END_TIME_COLNAME),
        (START_TIME_COLNAME, BG_START_TIME_COLNAME, BG_END_TIME_COLNAME, END_TIME_COLNAME),
    )

    ############### Particle Data ###############
    @property
    def DATA_PATH(self):
        return f"{getcwd()}/particle_data"

    PROTON_COLUMN_PREFIX = "H_Flux"
    ELECTRON_COLUMN_PREFIX = "Electron_Flux"
    PARTICLE_COLUMN_PREFIX = {"protons": PROTON_COLUMN_PREFIX, "electrons": ELECTRON_COLUMN_PREFIX}

    # Flux column names of the loaded data per sensor and particle
    RAW_FLUX_COLUMN = {
        "het": {"protons": "H_Flux", "electrons": "Electron_Flux"},
        "ept": {"protons": "Ion_Flux", "electrons": "Electron_Flux"},
    }

    # Energy bins keys of the loaded data per sensor and particle
    RAW_ENERGY_BINS_COLUMN = {
        "het": {"protons": "H_Bins", "electrons": "Electron_Bins"},
        "ept": {"protons": "Ion_Bins", "electrons": "Electron_Bins"},
    }

    ############### VDA ###############
    M_REST = {"protons": 938.27, "electrons": 0.511}

    def _epd_load(self, *args, **kwargs):
        return epd_load(*args, **kwargs)

    def _check_times_file_columns(self, filepath: str, columns: pd.Index) -> None:
        layouts = (self.REFERENCE_TIMES_FILE_LAYOUT,) + self.DATE_RANGE_FILE_LAYOUTS
        if set(columns) not in [set(layout) for layout in layouts]:
            accepted = "\n".join(f"    {self.EVENT_INDEX_NAME}, {', '.join(layout)}" for layout in layouts)
            raise ValueError(
                f"Unexpected columns in {filepath}: {', '.join(columns)}\n"
                f"Accepted columns:\n{accepted}"
            )

    def times_file_type(self, filepath: str) -> str:
        """Returns "reference times" or "datetime ranges", deduced from the columns of the events file"""
        columns = pd.read_csv(filepath, sep=",", header=0, index_col=0, skipinitialspace=True, nrows=0).columns.str.strip()
        self._check_times_file_columns(filepath, columns)
        return "reference times" if self.REF_TIME_COLNAME in columns else "datetime ranges"

    def _read_times_file(self, filepath: str) -> pd.DataFrame:
        """Reads the events file. Its type (datetime ranges or reference times) is deduced from its columns"""
        df = pd.read_csv(filepath, sep=",", header=0, index_col=0, skipinitialspace=True)
        df.columns = df.columns.str.strip()
        df.index.name = self.EVENT_INDEX_NAME
        self._check_times_file_columns(filepath, df.columns)
        df = df.apply(pd.to_datetime)
        if self.REF_TIME_COLNAME in df.columns:
            df[self.START_TIME_COLNAME] = df[self.REF_TIME_COLNAME] - timedelta(hours=self.parameters.bg_hours_prior)
            df[self.END_TIME_COLNAME] = df[self.REF_TIME_COLNAME] + timedelta(hours=self.parameters.bg_hours_after)
        elif self.BG_END_TIME_COLNAME in df.columns and self.BG_START_TIME_COLNAME not in df.columns:
            df[self.BG_START_TIME_COLNAME] = df[self.START_TIME_COLNAME]
        for col in (self.BG_START_TIME_COLNAME, self.BG_END_TIME_COLNAME, self.REF_TIME_COLNAME):
            if col not in df.columns:
                df[col] = pd.NaT
        return df

    def construct_times_df(self):
        """Creates self.df_times with the data range and background window of each event.

        The events are read from parameters.input_filepath, or if it is empty, from parameters.date_ranges.
        Events without a background window in the events file get the default one
        (parameters.bg_after_start minutes after the start time).
        """
        if self.parameters.input_filepath:
            self.df_times = self._read_times_file(self.parameters.input_filepath)
        else:
            date_ranges = self.parameters.date_ranges
            if len(date_ranges) == 0:
                raise ValueError("No events: set input_filepath or add a datetime range to date_ranges")
            for index_event, (start, end) in enumerate(date_ranges, start=1):
                if start is None or end is None or start >= end:
                    raise ValueError(f"Event {index_event}: the datetime range start ({start}) must be before its end ({end})")
            self.df_times = pd.DataFrame(
                {
                    self.START_TIME_COLNAME: [start for start, _ in date_ranges],
                    self.BG_START_TIME_COLNAME: [pd.NaT] * len(date_ranges),
                    self.BG_END_TIME_COLNAME: [pd.NaT] * len(date_ranges),
                    self.END_TIME_COLNAME: [end for _, end in date_ranges],
                    self.REF_TIME_COLNAME: [pd.NaT] * len(date_ranges),
                },
                index=pd.Index(range(1, len(date_ranges) + 1), name=self.EVENT_INDEX_NAME),
            )
        self.df_times = self.df_times[[
            self.START_TIME_COLNAME,
            self.BG_START_TIME_COLNAME,
            self.BG_END_TIME_COLNAME,
            self.END_TIME_COLNAME,
            self.REF_TIME_COLNAME,
        ]]

        for col in (self.BG_START_TIME_COLNAME, self.BG_END_TIME_COLNAME, self.REF_TIME_COLNAME):
            self.df_times[col] = self.df_times[col].astype(self.df_times[self.START_TIME_COLNAME].dtype)

        self._bg_sources = {
            index_event: "input file"
            if pd.notna(row[self.BG_START_TIME_COLNAME]) and pd.notna(row[self.BG_END_TIME_COLNAME])
            else "default"
            for index_event, row in self.df_times.iterrows()
        }
        self.set_default_bg_window(*self.parameters.bg_after_start)

        if self.parameters.view_dfs:
            return self.df_times

    def _bg_window(self, index_event) -> tuple:
        """Returns the (start, end) of the background window of the event"""
        event_times = self.df_times.loc[index_event]
        return (
            event_times[self.BG_START_TIME_COLNAME].to_pydatetime(),
            event_times[self.BG_END_TIME_COLNAME].to_pydatetime(),
        )

    def set_default_bg_window(self, start_minutes: int, end_minutes: int) -> None:
        """Sets the default background window, in minutes after the start time of the event.

        It applies to the events without a background window from the events file or set_bg_window.
        """
        self.parameters.bg_after_start = (start_minutes, end_minutes)
        for index_event, source in self._bg_sources.items():
            if source == "default":
                self._apply_default_bg_window(index_event)

    def _apply_default_bg_window(self, index_event) -> None:
        start_minutes, end_minutes = sorted(self.parameters.bg_after_start)
        start_time = self.df_times.loc[index_event, self.START_TIME_COLNAME]
        self.df_times.loc[index_event, self.BG_START_TIME_COLNAME] = start_time + timedelta(minutes=start_minutes)
        self.df_times.loc[index_event, self.BG_END_TIME_COLNAME] = start_time + timedelta(minutes=end_minutes)

    def reset_bg_window(self, index_event) -> None:
        """Sets the background window of the event back to the default one"""
        self._bg_sources[index_event] = "default"
        self._apply_default_bg_window(index_event)

    def _bg_window_points(self, index_event, times: pd.DatetimeIndex) -> int:
        bg_start, bg_end = self._bg_window(index_event)
        return int(((times >= bg_start) & (times <= bg_end)).sum())

    def _bg_window_warnings(self, index_event, times: pd.DatetimeIndex) -> list:
        """Warnings for a background window outside the data of the event or with less than 3 points"""
        bg_start, bg_end = self._bg_window(index_event)
        event_times = self.df_times.loc[index_event]
        warnings = []
        if bg_start < event_times[self.START_TIME_COLNAME] or bg_end > event_times[self.END_TIME_COLNAME]:
            warnings.append(
                f"Warning: the background window of event {index_event} ({bg_start} to {bg_end}) is outside "
                f"its data range ({event_times[self.START_TIME_COLNAME]} to {event_times[self.END_TIME_COLNAME]})"
            )
        n_points = self._bg_window_points(index_event, times)
        if n_points < 3:
            warnings.append(f"Warning: the background window of event {index_event} has {n_points} data points")
        return warnings

    def _check_bg_window(self, index_event, times: pd.DatetimeIndex) -> None:
        for warning in self._bg_window_warnings(index_event, times):
            print(warning)

    def _bg_window_source(self, index_event) -> str:
        return self._bg_sources[index_event]

    def set_bg_window(self, index_event, start, end) -> None:
        """Sets the background window of the event, used instead of the input file or default one"""
        start, end = pd.Timestamp(start), pd.Timestamp(end)
        if start >= end:
            raise ValueError(f"Background start ({start}) must be before its end ({end})")
        self.df_times.loc[index_event, self.BG_START_TIME_COLNAME] = start
        self.df_times.loc[index_event, self.BG_END_TIME_COLNAME] = end
        self._bg_sources[index_event] = "set"

    def save_times(self, filepath: str) -> None:
        """Saves the data range and background window of each event.

        The saved file can be used as parameters.input_filepath to reproduce the same windows.
        """
        self.df_times[[
            self.START_TIME_COLNAME,
            self.BG_START_TIME_COLNAME,
            self.BG_END_TIME_COLNAME,
            self.END_TIME_COLNAME,
        ]].to_csv(filepath, date_format="%Y-%m-%d %H:%M:%S")

    def construct_energies_df(self):
        """Creates self.df_energies with the energy bins of every channel of the available sensors and particles"""
        df_sensors = {}
        for sensor, particles in self.parameters.AVAILABLE_SENSORS_PARTICLES.items():
            if len(particles) == 0:
                continue
            df_protons, df_electrons, energies = self._epd_load(
                sensor=sensor,
                level="l2",
                startdate=self.df_times.iloc[0][self.START_TIME_COLNAME],
                enddate=self.df_times.iloc[0][self.END_TIME_COLNAME],
                viewing="sun",
                path=self.DATA_PATH,
                autodownload=True,
            )
            df_particles = []
            for particle, df_particle in (("protons", df_protons), ("electrons", df_electrons)):
                if particle not in particles:
                    continue
                particle_prefix = self.PARTICLE_COLUMN_PREFIX[particle]
                df_particle = df_particle.rename(
                    lambda x: x.replace(self.RAW_FLUX_COLUMN[sensor][particle], particle_prefix),
                    axis="columns",
                )
                energy_bins = self.RAW_ENERGY_BINS_COLUMN[sensor][particle]
                df_energies = pd.DataFrame(
                    {
                        "Low Energy": energies[f"{energy_bins}_Low_Energy"],
                        "Bin Width": energies[f"{energy_bins}_Width"],
                    },
                    index=df_particle[particle_prefix].columns,
                )
                df_energies["High Energy"] = df_energies["Low Energy"] + df_energies["Bin Width"]
                df_particles.append(df_energies)
            df_sensors[sensor] = pd.concat(df_particles)

        self.df_energies = pd.concat(df_sensors.values(), keys=list(df_sensors), names=["sensor", "channel"])

        if self.parameters.view_dfs:
            return self.df_energies

    def _iter_sensor_particles(self):
        """Yields (sensor, particle, particle_prefix) for the selected sensors and particles"""
        for sensor, particles in self.parameters.sensors_particles.items():
            for particle in particles:
                yield sensor, particle, self.PARTICLE_COLUMN_PREFIX[particle]

    def _iter_sensor_particle_viewings(self):
        """Yields (sensor, particle, viewing, particle_prefix) for the selected sensors, particles and viewings"""
        for sensor, particle, particle_prefix in self._iter_sensor_particles():
            for viewing in self.parameters.viewings:
                yield sensor, particle, viewing, particle_prefix

    def _download_data(self, show_progress: bool = True) -> pd.DataFrame:
        df_rows = []
        keys = []
        for index, row in self.df_times.iterrows():
            if show_progress:
                print(f"Working on event {index}...")
            df_row = pd.DataFrame({})
            keys.append(index)
            for sensor, particles in self.parameters.sensors_particles.items():
                
                if len(particles) == 0:
                    continue
                
                for viewing in self.parameters.viewings:
                    df_protons, df_electrons, _ = self._epd_load(
                        sensor=sensor,
                        level="l2",
                        startdate=row[self.START_TIME_COLNAME],
                        enddate=row[self.END_TIME_COLNAME],
                        viewing=viewing,
                        path=self.DATA_PATH,
                        autodownload=True,
                    )
                    for particle, df_particle in (("protons", df_protons), ("electrons", df_electrons)):
                        if particle not in particles:
                            continue
                        flux_cols_name = self.RAW_FLUX_COLUMN[sensor][particle]
                        df_particle = df_particle[
                            [c for c in df_particle.columns if c[0] == flux_cols_name]
                        ]
                        df_particle = df_particle[
                            (df_particle.index >= row[self.START_TIME_COLNAME])
                            & (df_particle.index <= row[self.END_TIME_COLNAME])
                        ]
                        if self.parameters.resample_frequency:
                            df_particle = df_particle.resample(
                                self.parameters.resample_frequency, origin="start"
                            ).mean()
                            df_particle.index = df_particle.index.floor("min")
                        df_particle = pd.concat(
                            [df_particle],
                            keys=[(sensor, particle, viewing)],
                            axis="columns",
                        )
                        df_particle = df_particle.rename(
                            lambda x: x.replace(
                                flux_cols_name, self.PARTICLE_COLUMN_PREFIX[particle]
                            ),
                            axis="columns",
                        )
                        df_row = pd.concat([df_row, df_particle], axis="columns")
            df_rows.append(df_row)

        if show_progress:
            print(f"Done")
        return pd.concat(df_rows, keys=keys, names=[self.EVENT_INDEX_NAME, "Time"])

    def construct_particles_df(self):
        if self.parameters.load_data:
            self.df_data = pd.read_pickle(self.parameters.load_data_filepath)
        else:
            self.df_data = self._download_data()
            if self.parameters.save_data:
                self.df_data.to_pickle(self.parameters.save_data_filepath)

        if self.parameters.view_dfs:
            return self.df_data

    def _group_channels_de(
        self,
        df_to_group: pd.DataFrame,
        groups: list[list[str]],
        energy_bins_widths: list[list[float]],
        names: list[str]
    ) -> pd.DataFrame:
        # I = ΣI_n*ΔE_n / ΣΔE_n
        
        grouped_all = {}
        for columns, energy_bins_width, name in zip(groups, energy_bins_widths, names):
            de = sum(energy_bins_width)
            grouped_series = df_to_group.loc[slice(None), columns[0]]*energy_bins_width[0]
            for column, eb_width in zip(columns[1:], energy_bins_width[1:]):
                grouped_series = grouped_series.add(df_to_group.loc[slice(None), column]*eb_width, fill_value=0)
            grouped_all[name] = grouped_series/de
        df_grouped = pd.DataFrame(grouped_all)
        return df_grouped

    def group_energy_channels(self):
        grouped_frames = []
        for sensor, particle, viewing, particle_prefix in self._iter_sensor_particle_viewings():
            df_temp = self._group_channels_de(
                self.df_data[sensor][particle][viewing][particle_prefix],
                [[f"{particle_prefix}_{c}" for c in spec["channels"]]
                 for spec in self.parameters.channel_groups[particle].values()
                 if spec["sensor"] == sensor],
                [[self.df_energies.loc[(sensor, f"{particle_prefix}_{c}"), "Bin Width"] for c in spec["channels"]]
                 for spec in self.parameters.channel_groups[particle].values()
                 if spec["sensor"] == sensor],
                [key for key, spec in self.parameters.channel_groups[particle].items() if spec["sensor"] == sensor]
            )
            grouped_frames.append(pd.concat(
                [df_temp],
                keys=[(sensor, particle, viewing, particle_prefix)],
                axis="columns",
            ))
        self.df_grouped = pd.concat(grouped_frames, axis="columns")

        if self.parameters.view_dfs:
            return self.df_grouped

    def _onset_detection_sigma(
        self,
        series: pd.Series,
        s: int = 3,
        n: int = 3,
        bg_start: int | datetime = 0,
        bg_end: int | datetime = 12,
    ) -> tuple:
        """Returns:

        1. Onset time or None if no event detected
        2. Background start
        3. Background end
        4. Background Level
        5. Threshold
        """
        if isinstance(bg_start, (int, np.integer)):
            bg_start = series.index[bg_start]
        if isinstance(bg_end, (int, np.integer)):
            bg_end = series.index[bg_end]
        bg_level = (bg_series := series[bg_start:bg_end]).mean()
        threshold = bg_level + s * bg_series.std()
        onset_time = None

        streak = 0
        for index, value in series.items():
            if value > threshold:
                streak += 1
                if onset_time is None:
                    onset_time = index
            else:
                streak = 0
                onset_time = None

            if streak >= n:
                break

        if streak < n:
            onset_time = None

        return (
            onset_time,
            bg_start,
            bg_end,
            {"bg_level": bg_level, "threshold": threshold},
        )

    def _onset_detection_poisson_cusum_bootstrap(
        self,
        series: pd.Series,
        sensor: str,
        particle: str,
        viewing: str,
        channel: str,
        bg_start: int | datetime = 0,
        bg_end: int | datetime = 12,
        bootstraps: int = 1000,
        cusum_minutes: int = 60,
        sample_size: float = 0.75,
        limit_averaging: str = "4 min",
    ) -> tuple:
        if isinstance(bg_start, (int, np.integer)):
            bg_start = series.index[bg_start]
        if isinstance(bg_end, (int, np.integer)):
            bg_end = series.index[bg_end]
        df = pd.DataFrame(series)
        df.index.freq = self.parameters.resample_frequency
        protons = Onset(
            spacecraft="Solar Orbiter",
            sensor=sensor.upper(),
            species=particle[
                :-1
            ],  # rest of the tool uses the plurar form. This function needs singular, so omit the final "s"
            viewing=viewing,
            data_level="l2",
            data_path="",
            start_date="",
            end_date="",
            data=df,
        )
        channels = channel.split("-")
        protons.set_custom_channel_energies(
            low_bounds=[self.df_energies.loc[sensor, channels[0]]["Low Energy"]],
            high_bounds=[self.df_energies.loc[sensor, channels[1]]["High Energy"]],
            unit="MeV",
        )
        bg = BootstrapWindow(
            start=bg_start.strftime("%Y-%m-%d %H:%M"),
            end=bg_end.strftime("%Y-%m-%d %H:%M"),
            bootstraps=bootstraps,
        )
        rng = 101010101
        protons.onset_statistics_per_channel(
            channels=channel,
            background=bg,
            cusum_minutes=cusum_minutes,
            sample_size=sample_size,
            viewing=protons.viewing,
            limit_averaging=limit_averaging,
            random_seed=rng,
            print_output=False,
        )

        return (
            protons.onset_statistics[channel][0],
            bg_start,
            bg_end,
            protons.onset_statistics[channel],
        )

    def _onset_detection(
        self, series: pd.Series, method: str = "sigma", **kwargs
    ) -> tuple:
        if method == "sigma":
            onset_results = self._onset_detection_sigma(
                series, kwargs["s"], kwargs["n"], kwargs["bg_start"], kwargs["bg_end"]
            )
        elif method == "poisson_cusum_bootstrap":
            onset_results = self._onset_detection_poisson_cusum_bootstrap(
                series,
                kwargs["sensor"],
                kwargs["particle"],
                kwargs["viewing"],
                kwargs["channel"],
                kwargs["bg_start"],
                kwargs["bg_end"],
                kwargs["bootstraps"],
                kwargs["cusum_minutes"],
                kwargs["sample_size"],
                kwargs["limit_averaging"],
            )
        else:
            raise ValueError(f'Method named "{method}" is not implented')
        return onset_results

    def _onset_detection_df(
        self, df: pd.DataFrame, method: str = "sigma", **kwargs
    ) -> pd.DataFrame:
        rows = []
        index = []
        for index_event, df_event in df.groupby(level=0):
            self._check_bg_window(index_event, df_event.index.droplevel(0))
            bg_start, bg_end = self._bg_window(index_event)
            for sensor, particle, viewing, particle_prefix in self._iter_sensor_particle_viewings():
                df_inner = df_event[sensor][particle][viewing][particle_prefix]
                for column_name in df_inner.columns:
                    new_kwargs = dict(kwargs)
                    new_kwargs["sensor"] = sensor
                    new_kwargs["particle"] = particle
                    new_kwargs["viewing"] = viewing
                    new_kwargs["channel"] = column_name
                    new_kwargs["bg_start"] = bg_start
                    new_kwargs["bg_end"] = bg_end
                    try:
                        onset_time, bg_start, bg_stop, method_specific = self._onset_detection(
                            df_inner[column_name].droplevel(0, axis="index"),
                            method,
                            **new_kwargs,
                        )
                    except Exception as e:
                        print(index_event, type(e).__name__, new_kwargs)
                        onset_time, bg_start, bg_stop, method_specific = pd.NaT, pd.NaT, pd.NaT, None
                    rows.append({
                        "Onset Time": onset_time,
                        "Background Start": bg_start,
                        "Background End": bg_stop,
                        "Method Specific": method_specific,
                    })
                    index.append((index_event, sensor, particle, viewing, particle_prefix, column_name))

        return pd.DataFrame(
            rows,
            index=pd.MultiIndex.from_tuples(
                index,
                names=[self.EVENT_INDEX_NAME, "sensor", "particle", "viewing", "prefix", "channels"],
            ),
        )

    def calculate_onsets(self):
        for name in ("bg_start", "bg_end"):
            if name in self.parameters.onset_method_parameters:
                raise ValueError(
                    f"The {name} onset method parameter was removed in v0.3.0. The default background window "
                    f"is set with set_default_bg_window (minutes after the start time), "
                    f"and the window of an event with set_bg_window."
                )
        self.df_onsets = self._onset_detection_df(
            self.df_grouped,
            self.parameters.onset_method,
            **self.parameters.onset_method_parameters,
        )

        if self.parameters.view_dfs:
            return self.df_onsets

    def clean_onsets(self):
        self.df_onsets_existing = self.df_onsets[~pd.isna(self.df_onsets["Onset Time"])]

        if self.parameters.view_dfs:
            return self.df_onsets_existing

    def construct_options_df(self):
        self.df_options = self.df_onsets_existing.reorder_levels(
            [
                self.EVENT_INDEX_NAME,
                "sensor",
                "particle",
                "prefix",
                "channels",
                "viewing",
            ]
        )

        if self.parameters.view_dfs:
            return self.df_options

    def _channel_energy_range(self, sensor, particle, particle_prefix, channel) -> tuple:
        """Returns the (low, high) energy of a grouped channel in MeV"""
        channels = self.parameters.channel_groups[particle][channel]["channels"]
        low_energy = self.df_energies.loc[(sensor, f"{particle_prefix}_{channels[0]}"), "Low Energy"]
        high_energy = self.df_energies.loc[(sensor, f"{particle_prefix}_{channels[-1]}"), "High Energy"]
        return low_energy, high_energy

    def construct_energy_channels_characteristics(self):
        rows = []
        index = []
        for sensor, particle, particle_prefix in self._iter_sensor_particles():
            for channel in list(
                self.df_grouped[sensor][particle][self.parameters.viewings[0]][
                    particle_prefix
                ].columns
            ):
                low_energy, high_energy = self._channel_energy_range(sensor, particle, particle_prefix, channel)
                geo_mean = sqrt(low_energy) * sqrt(high_energy)
                inv_beta = 1 / sqrt(
                    1 - (1 / (1 + geo_mean / self.M_REST[particle])) ** 2
                )
                rows.append({"Geomagnetic Mean": geo_mean, "Inverse Beta": inv_beta})
                index.append((sensor, particle, channel))
        self.df_channels_chars = pd.DataFrame(
            rows,
            index=pd.MultiIndex.from_tuples(index, names=["sensor", "particle", "channel"]),
        )

        if self.parameters.view_dfs:
            return self.df_channels_chars

    def define_spacecraft_parameters(self):
        kernel_urls = [
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
        ]
        kernel_urls = [f"https://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/{url}"
                    for url in kernel_urls]

        kernel_files = [cache.download(url) for url in kernel_urls]

        spice.initialize(kernel_files)

    @staticmethod
    def _format_timedelta(td) -> str:
        return str(pd.Timedelta(td).to_pytimedelta()).split(".")[0]

    def print_results(self, events=None) -> None:
        """Prints the VDA results of the given event number(s), or of all events if None."""
        if events is None:
            events = self.results.index
        elif not isinstance(events, (list, tuple, pd.Index)):
            events = [events]

        for index_event in events:
            time_start = self.df_times.loc[index_event][self.START_TIME_COLNAME].strftime("%Y-%m-%d %H:%M")
            time_end = self.df_times.loc[index_event][self.END_TIME_COLNAME].strftime("%Y-%m-%d %H:%M")
            print(f"Event {index_event} ({time_start} to {time_end})")
            res = self.results.loc[index_event]
            if pd.isna(res["APL"]):
                print("    No results (not enough onset points)\n")
                continue
            print(f"    Release Time : {res['Release Time']} ± {self._format_timedelta(res['Release Time Error'])}")
            print(f"    Extra Time   : {self._format_timedelta(res['Extra Time'])}")
            print(f"    APL          : {res['APL']:.2f} ± {res['APL Error']:.2f}\n")

    def compute_vda(self):
        """Fits the VDA line of each event and stores the results in self.results"""
        self._vda_fits = {}
        for index_event in self.df_options.index.unique(level=0):
            vda_points = []
            # onset times are fitted in seconds from the event start
            t0 = self.df_times.loc[index_event][self.START_TIME_COLNAME].to_pydatetime()
            t_sun_to_observer = (
                spice.get_body(
                    "Solar Orbiter",
                    self.df_times.loc[index_event][self.START_TIME_COLNAME],
                    spice_frame="SOLO_HEEQ"
                )
                .distance
                / const.c
            ).to(u.s).value
            for i, row in self.parameters.selected_onsets.loc[index_event].iterrows():
                if row["Viewing"] is None:
                    continue
                sensor, particle, particle_prefix, channel = i
                vda_points.append(
                    (self.df_channels_chars.loc[
                        sensor, particle, channel
                    ]["Inverse Beta"],
                    (self.df_onsets_existing.loc[
                        index_event,
                        sensor,
                        particle,
                        row["Viewing"],
                        particle_prefix,
                        channel
                    ]["Onset Time"]
                    .to_pydatetime() - t0)
                    .total_seconds())
                )

            if len(vda_points) < 2:
                # Not enough points for the linear regression
                # Consider throughing warning
                print(f"Not enough onset points in event {index_event}.")
                self.results.loc[index_event] = np.nan
                continue

            vda_points = sorted(vda_points, key=lambda x: x[0])
            inv_betas = np.array([p[0] for p in vda_points])
            onset_seconds = np.array([p[1] for p in vda_points])

            try:
                p, V = np.polyfit(inv_betas, onset_seconds, 1, cov=True)
                a = p[0]
                b = p[1]
                a_error = np.sqrt(V[0][0])
                b_error = np.sqrt(V[1][1])
            except ValueError:
                # not enough points for cov matrix
                print(
                    f"Not enough points for covariance matrix generation in event {index_event}"
                )
                a, b = np.polyfit(inv_betas, onset_seconds, 1)
                a_error = 0
                b_error = 0

            self.results.loc[index_event] = {
                "Release Time": (t0 + timedelta(seconds=b + t_sun_to_observer)).strftime('%Y-%m-%d %H:%M:%S'),
                "Release Time Error": timedelta(seconds=b_error),
                "Extra Time": timedelta(seconds=t_sun_to_observer),
                "APL": a / t_sun_to_observer,
                "APL Error": a_error / t_sun_to_observer,
            }
            self._vda_fits[index_event] = {
                "t0": t0,
                "inv_betas": inv_betas,
                "onset_seconds": onset_seconds,
                "a": a,
                "b": b,
                "b_error": b_error,
            }

    def plot_vda(self, savefig: bool = True, returnfig: bool = False):
        """Prints the results and plots the VDA fit of each event computed by compute_vda.

        If returnfig is True, returns the figure, or a list of figures if more than one event was plotted.
        """
        figs = []
        for index_event, fit in self._vda_fits.items():
            inv_betas = fit["inv_betas"]
            a, b, b_error = fit["a"], fit["b"], fit["b_error"]

            def to_time(seconds):
                return fit["t0"] + timedelta(seconds=seconds)

            res = self.results.loc[index_event]
            self.print_results(index_event)

            fig, ax = plt.subplots(figsize=(12, 8), layout="constrained")
            ax.scatter(
                inv_betas,
                [to_time(t) for t in fit["onset_seconds"]],
                color="black",
            )
            ax.plot(
                inv_betas,
                [to_time(a * x + b) for x in inv_betas],
                label="Linear Regression",
                color="blue",
            )
            ax.fill_between(
                inv_betas,
                [to_time(a * x + b - 2 * b_error) for x in inv_betas],
                [to_time(a * x + b + 2 * b_error) for x in inv_betas],
                color="blue",
                alpha=0.1,
            )
            fig.suptitle(f"Event {index_event} ({self.df_grouped.loc[index_event].index[1].to_pydatetime().strftime('%Y-%m-%d')})")
            ax.set_xlabel("Inverse Beta")
            ax.set_ylabel("Time")
            time_formatter = mdates.DateFormatter("%H:%M")
            ax.yaxis.set_major_formatter(time_formatter)
            ax.plot(
                [],
                [],
                alpha=0,
                label=f"Extra Time = {self._format_timedelta(res['Extra Time'])}",
            )
            ax.plot(
                [],
                [],
                alpha=0,
                label=f"Release Time = {res['Release Time']} +/- {self._format_timedelta(res['Release Time Error'])}",
            )
            ax.plot(
                [],
                [],
                alpha=0,
                label=f"APL = {res['APL']:.2f} +/- {res['APL Error']:.2f}",
            )
            # legend between the title and the plot
            ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=2, frameon=False)
            if savefig:
                # date_str = self.df_grouped.loc[index_event].index[1].to_pydatetime().strftime('%Y-%m-%d')
                time_start_str = self.df_times.loc[index_event][self.START_TIME_COLNAME].strftime("%Y-%m-%d_%H%M")
                time_end_str = self.df_times.loc[index_event][self.END_TIME_COLNAME].strftime("%Y-%m-%d_%H%M")
                date_str = f"{time_start_str}_{time_end_str}"
                particles_str = "_".join([f'{s}-{p}' for s, ps in self.parameters.sensors_particles.items() for p in ps])
                freq_str = self.parameters.resample_frequency if self.parameters.resample_frequency != "" else "noresample"
                filename = f"{date_str}_{particles_str}_{freq_str}.png"
                plt.savefig(filename)
            plt.show()
            figs.append(fig)

        if returnfig:
            return figs[0] if len(figs) == 1 else figs

    def plot(self, savefig: bool = True, returnfig: bool = False):
        """Fits the VDA line of each event, stores it in self.results, prints the results and plots them.

        If returnfig is True, returns the figure, or a list of figures if more than one event was plotted.
        """
        self.compute_vda()
        return self.plot_vda(savefig, returnfig)

    def _plot_event_bg(self, ax, event_no) -> None:
        """Plots the grouped channels of the event with its background window"""
        # one colormap per sensor and particle, darker for the later channels of the group
        colormaps = ["Blues", "Oranges", "Greens", "Purples", "Reds", "Greys"]
        linestyles = ["-", "--", ":", "-."]
        viewings = self.parameters.viewings
        # zeros cannot be shown in log scale
        temp_df = self.df_grouped.loc[event_no].replace(0, np.nan)
        for i_group, (sensor, particle, particle_prefix) in enumerate(self._iter_sensor_particles()):
            cmap = plt.get_cmap(colormaps[i_group % len(colormaps)])
            for i_viewing, viewing in enumerate(viewings):
                df_channels = temp_df[sensor][particle][viewing][particle_prefix]
                n_channels = len(df_channels.columns)
                for i_channel, channel in enumerate(df_channels.columns):
                    low_energy, high_energy = self._channel_energy_range(sensor, particle, particle_prefix, channel)
                    label = f"{sensor.upper()} {particle} {low_energy:.2f}-{high_energy:.2f} MeV"
                    if len(viewings) > 1:
                        label += f" ({viewing})"
                    ax.plot(
                        df_channels[channel],
                        color=cmap(0.4 + 0.6 * i_channel / max(n_channels - 1, 1)),
                        linestyle=linestyles[i_viewing % len(linestyles)],
                        label=label,
                    )

        bg_start, bg_end = self._bg_window(event_no)
        ax.axvspan(bg_start, bg_end, color="green", alpha=0.2, label="Background")
        bg_end_format = "%H:%M" if bg_start.date() == bg_end.date() else "%Y-%m-%d %H:%M"
        ax.set_title(
            f"Event {event_no}\n"
            f"Background: {bg_start:%Y-%m-%d %H:%M} to {bg_end:{bg_end_format}} ({self._bg_window_source(event_no)})"
        )

        ax.set_yscale("log")
        ax.set_ylabel(r"Intensity (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$)")
        ax.set_xlim(temp_df.index[0], temp_df.index[-1])
        locator = mdates.AutoDateLocator()
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        ax.set_xlabel("Time (UTC)")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize="small")

    def _plot_channel_onsets(self, event_no, sensor, particle, particle_prefix, channel, selected_viewing=None):
        """Plots the detected onsets of a grouped channel of the event, one subplot per viewing"""
        viewings = self.parameters.viewings
        temp_df = self.df_grouped.loc[event_no]
        ncols = min(len(viewings), 3)
        nrows = ceil(len(viewings) / ncols)
        fig = plt.figure(figsize=(14, 4.5 * nrows + 1), dpi=300, layout="constrained")
        # the last row holds the legend of all the viewings
        grid = fig.add_gridspec(nrows + 1, ncols, height_ratios=[1] * nrows + [0.15])
        axs_flat = [fig.add_subplot(grid[row, col]) for row in range(nrows) for col in range(ncols)]
        ax_legend = fig.add_subplot(grid[nrows, :])
        ax_legend.axis("off")
        for ax in axs_flat[len(viewings):]:
            ax.axis("off")
        for ax, viewing in zip(axs_flat, viewings):
            ax.set_title(f"{viewing} (selected)" if viewing == selected_viewing else viewing,
                         fontweight="bold" if viewing == selected_viewing else "normal")
            try:
                onset_results = self.df_onsets_existing.loc[(event_no, sensor, particle, viewing, particle_prefix, channel)]
            except KeyError:
                ax.text(0.5, 0.5, "No onset", transform=ax.transAxes, ha="center", va="center")
                continue

            ax.plot(temp_df[sensor][particle][viewing][particle_prefix][channel].fillna(0).ffill(), label="Data")
            ax.set_yscale("log")
            ax.axvspan(onset_results["Background Start"], onset_results["Background End"],
                       color="green", alpha=0.3, label="BG sample")
            values = [f"Onset {onset_results['Onset Time']:%H:%M}"]
            # bg level and threshold are only provided by the sigma method
            method_specific = onset_results["Method Specific"]
            if isinstance(method_specific, dict) and "bg_level" in method_specific:
                ax.axhline(method_specific["bg_level"], color="green", linestyle="dashed", label="BG level")
                ax.axhline(method_specific["threshold"], color="red", linestyle="dashed", label="Threshold")
                values += [f"BG {method_specific['bg_level']:.3g}", f"Threshold {method_specific['threshold']:.3g}"]
            ax.axvline(onset_results["Onset Time"], color="purple", linestyle="dashed", label="Onset")
            ax.text(0.98, 0.03, "\n".join(values), transform=ax.transAxes, ha="right", va="bottom",
                    fontsize="x-small", bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})

            ax.set_xlim(temp_df.index[0], temp_df.index[-1])
            locator = mdates.AutoDateLocator(maxticks=5)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))

        # one legend for all the viewings, from the first subplot with an onset
        for ax in axs_flat:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax_legend.legend(handles, labels, loc="center", ncols=len(handles), fontsize="small", frameon=False)
                break

        low_energy, high_energy = self._channel_energy_range(sensor, particle, particle_prefix, channel)
        fig.suptitle(
            f"Detected onsets for event {event_no} ({temp_df.index[0]:%Y-%m-%d}) | "
            f"{sensor}/{particle} ({low_energy:.2f}-{high_energy:.2f} MeV)"
        )
        return fig

    def plot_bg_selection(self):
        """Plots the grouped channels of each event with its background window"""
        for event_no in self.df_grouped.index.unique(level=0):
            self._check_bg_window(event_no, self.df_grouped.loc[event_no].index)
            fig, ax = plt.subplots(figsize=(12, 6))
            self._plot_event_bg(ax, event_no)
            fig.tight_layout()
            plt.show()
