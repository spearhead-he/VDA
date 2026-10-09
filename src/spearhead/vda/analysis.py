"""Velocity Dispersion Analysis (VDA) of Solar Energetic Particle events.

The VDA class performs the steps of the analysis:
- events: data ranges and background windows, from an events file or datetime ranges
- data: download (or load) of the intensities of the observer, resampling and grouping of the energy channels
- onsets: onset times of each grouped channel and viewing, and the selection of the ones used for the fit
- VDA fit: release time and apparent path length of each event

Its parameters are in conf.py, the observers (spacecraft) in observers.py and its plots in views.py.
"""
import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const

from contextlib import redirect_stdout
from io import StringIO
from math import sqrt
from os import getcwd
from datetime import datetime, timedelta

from pyonset import onset_determination

from .observers import OBSERVERS


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

    ############### VDA ###############
    M_REST = {"protons": 938.27, "electrons": 0.511}

    @property
    def observer(self):
        return OBSERVERS[self.parameters.observer]

    def _channel_column(self, particle, channel) -> str:
        """Column name of the channel number in the data, e.g. H_Flux_12"""
        return f"{self.PARTICLE_COLUMN_PREFIX[particle]}_{channel}"

    def _load_observer_data(self, sensor, startdate, enddate, viewing, particles) -> dict:
        """Data of the particles of the sensor, with the channel column names of the data"""
        data = self.observer.load(sensor, startdate, enddate, viewing, self.DATA_PATH, particles=list(particles))
        return {
            particle: (
                particle_data.flux.rename(columns=lambda c: self._channel_column(particle, c)),
                particle_data.energies.rename(index=lambda c: self._channel_column(particle, c)),
            )
            for particle, particle_data in data.items()
        }

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

    def bg_window(self, index_event) -> tuple:
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
        # before construct_times_df there are no events yet, the default is applied when they are created
        for index_event, source in getattr(self, "_bg_sources", {}).items():
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

    def bg_window_points(self, index_event, times: pd.DatetimeIndex) -> int:
        bg_start, bg_end = self.bg_window(index_event)
        return int(((times >= bg_start) & (times <= bg_end)).sum())

    def bg_window_warnings(self, index_event, times: pd.DatetimeIndex) -> list:
        """Warnings for a background window outside the data of the event or with less than 3 points"""
        bg_start, bg_end = self.bg_window(index_event)
        event_times = self.df_times.loc[index_event]
        warnings = []
        if bg_start < event_times[self.START_TIME_COLNAME] or bg_end > event_times[self.END_TIME_COLNAME]:
            warnings.append(
                f"Warning: the background window of event {index_event} ({bg_start} to {bg_end}) is outside "
                f"its data range ({event_times[self.START_TIME_COLNAME]} to {event_times[self.END_TIME_COLNAME]})"
            )
        n_points = self.bg_window_points(index_event, times)
        if n_points < 3:
            warnings.append(f"Warning: the background window of event {index_event} has {n_points} data points")
        return warnings

    def check_bg_window(self, index_event, times: pd.DatetimeIndex) -> None:
        for warning in self.bg_window_warnings(index_event, times):
            print(warning)

    def bg_window_source(self, index_event) -> str:
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

    def construct_energies_df(self, startdate=None, enddate=None):
        """Creates self.df_energies with the energy bins of every channel of the available sensors and particles.

        The energy bins are read from the data of startdate to enddate, by default the data range of the first event.
        """
        if startdate is None or enddate is None:
            startdate = self.df_times.iloc[0][self.START_TIME_COLNAME]
            enddate = self.df_times.iloc[0][self.END_TIME_COLNAME]
        df_sensors = {}
        for sensor, particles in self.parameters.AVAILABLE_SENSORS_PARTICLES.items():
            if len(particles) == 0:
                continue
            # the energy bins are the same in all viewings
            data = self._load_observer_data(sensor, startdate, enddate, self.observer.SENSOR_VIEWINGS[sensor][0], particles)
            df_particles = []
            for particle, (_, df_energies) in data.items():
                df_energies = df_energies.copy()
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

    def sensor_viewings(self, sensor) -> list:
        """Selected viewings of the sensor, in the order of parameters.viewings"""
        return [v for v in self.parameters.viewings if v in self.observer.SENSOR_VIEWINGS[sensor]]

    def check_viewings(self) -> None:
        """Raises an error if a sensor of the grouped channels has none of the selected viewings"""
        for sensor in self.parameters.sensors_particles:
            if not self.sensor_viewings(sensor):
                raise ValueError(
                    f"None of the selected viewings ({', '.join(self.parameters.viewings)}) is a viewing of "
                    f"{self.observer.label} {sensor.upper()}: select one of "
                    f"{', '.join(self.observer.SENSOR_VIEWINGS[sensor])}, or remove its grouped channels"
                )

    def _iter_sensor_particle_viewings(self):
        """Yields (sensor, particle, viewing, particle_prefix) for the selected sensors, particles and viewings"""
        for sensor, particle, particle_prefix in self._iter_sensor_particles():
            for viewing in self.sensor_viewings(sensor):
                yield sensor, particle, viewing, particle_prefix

    def _download_data(self, show_progress: bool = True) -> pd.DataFrame:
        self.check_viewings()
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
                
                for viewing in self.sensor_viewings(sensor):
                    data = self._load_observer_data(
                        sensor, row[self.START_TIME_COLNAME], row[self.END_TIME_COLNAME], viewing, particles
                    )
                    for particle, (df_particle, _) in data.items():
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
                            keys=[(sensor, particle, viewing, self.PARTICLE_COLUMN_PREFIX[particle])],
                            axis="columns",
                        )
                        df_row = pd.concat([df_row, df_particle], axis="columns")
            df_rows.append(df_row)

        if show_progress:
            print(f"Done")
        return pd.concat(df_rows, keys=keys, names=[self.EVENT_INDEX_NAME, "Time"])

    def _use_saved_resample_frequency(self) -> None:
        def describe(frequency):
            return repr(frequency) if frequency else "no resampling"

        saved = self.df_data.attrs.get("resample_frequency")
        current = self.parameters.resample_frequency
        if saved is None:
            print(f"The loaded data do not include their resample frequency, the current one ({describe(current)}) is used")
        elif saved != current:
            self.parameters.resample_frequency = saved
            print(f"Resample frequency set to {describe(saved)}, the one of the loaded data (instead of {describe(current)})")

    def _check_saved_observer(self) -> None:
        # data saved before the observer was recorded are Solar Orbiter data
        saved = self.df_data.attrs.get("observer", "solo")
        if saved != self.parameters.observer:
            raise ValueError(
                f"The data of {self.parameters.load_data_filepath} are data of the {saved!r} observer, "
                f"but parameters.observer is {self.parameters.observer!r}"
            )

    def construct_particles_df(self):
        """Loads the data from parameters.load_data_filepath, or downloads them if it is empty.

        The data are saved to parameters.save_data_filepath, if it is set.
        """
        if self.parameters.load_data_filepath:
            self.df_data = pd.read_pickle(self.parameters.load_data_filepath)
            self._check_saved_observer()
            self._use_saved_resample_frequency()
        else:
            self.df_data = self._download_data()
            # saved with the data, so that loaded data are used with their observer and resample frequency
            self.df_data.attrs["observer"] = self.parameters.observer
            self.df_data.attrs["resample_frequency"] = self.parameters.resample_frequency
        if self.parameters.save_data_filepath:
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
        self.check_viewings()
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

    @staticmethod
    def _bg_times(series: pd.Series, bg_start, bg_end) -> tuple:
        """The background start and end as times, also when given as point indices"""
        if isinstance(bg_start, (int, np.integer)):
            bg_start = series.index[bg_start]
        if isinstance(bg_end, (int, np.integer)):
            bg_end = series.index[bg_end]
        return bg_start, bg_end

    def _onset_detection_sigma(
        self,
        series: pd.Series,
        s: int = 3,
        n: int = 3,
        bg_start: int | datetime = 0,
        bg_end: int | datetime = 12,
    ) -> tuple:
        """The onset is the last point before the first of n consecutive points above the threshold
        (as the Poisson-CUSUM method of pyonset).

        Returns:

        1. Onset time or None if no event detected
        2. Background start
        3. Background end
        4. Method specific values: background level and threshold
        """
        bg_start, bg_end = self._bg_times(series, bg_start, bg_end)
        bg_level = (bg_series := series[bg_start:bg_end]).mean()
        threshold = bg_level + s * bg_series.std()
        onset_time = None

        streak = 0
        previous = None
        for index, value in series.items():
            if value > threshold:
                streak += 1
                if onset_time is None:
                    # the first point of the series has no point before it
                    onset_time = index if previous is None else previous
            else:
                streak = 0
                onset_time = None
            previous = index

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

    def _onset_detection_poisson_cusum(
        self,
        series: pd.Series,
        cusum_minutes: int = 30,
        sigma_multiplier: int = 2,
        bg_start: int | datetime = 0,
        bg_end: int | datetime = 12,
    ) -> tuple:
        """Poisson-CUSUM onset determination of pyonset (onset_determination), after the background window.

        Returns:

        1. Onset time or None if no event detected
        2. Background start
        3. Background end
        4. Method specific values: background level, mu_d (background level + sigma_multiplier * standard
           deviation), and the k and h parameters of the CUSUM function
        """
        bg_start, bg_end = self._bg_times(series, bg_start, bg_end)
        bg_series = series[bg_start:bg_end]
        # the CUSUM window in data points, from the cadence of the series (also without resampling)
        cadence = series.index.to_series().diff().median()
        cusum_window = max(1, round(pd.Timedelta(minutes=cusum_minutes) / cadence))
        # pyonset prints its warnings (e.g. for a background of zeros)
        with redirect_stdout(StringIO()):
            bg_level, mu_d, k, h, _, _, onset_time = onset_determination(
                (np.nanmean(bg_series), np.nanstd(bg_series)),
                series,
                cusum_window,
                bg_end,
                sigma_multiplier=sigma_multiplier,
            )

        return (
            None if pd.isna(onset_time) else onset_time,
            bg_start,
            bg_end,
            {"bg_level": bg_level, "mu_d": mu_d, "k": k, "h": h},
        )

    # Onset determination method of each name of AVAILABLE_ONSET_METHODS
    ONSET_DETECTION_METHODS = {
        "sigma": "_onset_detection_sigma",
        "poisson_cusum": "_onset_detection_poisson_cusum",
    }

    def _onset_detection(self, series: pd.Series, method: str = "sigma", **kwargs) -> tuple:
        try:
            detect = getattr(self, self.ONSET_DETECTION_METHODS[method])
        except KeyError:
            raise ValueError(f"Unknown onset method {method!r}") from None
        return detect(series, **kwargs)

    def _onset_detection_df(
        self, df: pd.DataFrame, method: str = "sigma", **kwargs
    ) -> pd.DataFrame:
        rows = []
        index = []
        for index_event, df_event in df.groupby(level=0):
            self.check_bg_window(index_event, df_event.index.droplevel(0))
            bg_start, bg_end = self.bg_window(index_event)
            for sensor, particle, viewing, particle_prefix in self._iter_sensor_particle_viewings():
                df_inner = df_event[sensor][particle][viewing][particle_prefix]
                for column_name in df_inner.columns:
                    new_kwargs = dict(kwargs)
                    new_kwargs["bg_start"] = bg_start
                    new_kwargs["bg_end"] = bg_end
                    # the background window of the event is kept for the next channels, also after a failure
                    try:
                        onset_time, used_bg_start, used_bg_end, method_specific = self._onset_detection(
                            df_inner[column_name].droplevel(0, axis="index"),
                            method,
                            **new_kwargs,
                        )
                    except Exception as e:
                        print(index_event, type(e).__name__, new_kwargs)
                        onset_time, used_bg_start, used_bg_end, method_specific = pd.NaT, pd.NaT, pd.NaT, None
                    rows.append({
                        "Onset Time": onset_time,
                        "Background Start": used_bg_start,
                        "Background End": used_bg_end,
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
        method = self.parameters.onset_method
        expected = set(self.parameters.AVAILABLE_ONSET_METHODS[method])
        given = set(self.parameters.onset_method_parameters)
        if given - expected:
            raise ValueError(
                f"Parameters {sorted(given - expected)} are not parameters of the {method} onset method, "
                f"whose parameters are {sorted(expected)}"
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

    def select_onsets(self):
        """Selects the viewing of the onset used for each grouped channel ("Use all").

        The viewing is the first one with an onset, in the order of parameters.viewings. The interactive
        selection starts from the same viewings.
        """
        temp_df = self.df_options.droplevel(level=5)
        df_index = temp_df.index[~temp_df.index.duplicated(keep="first")]
        self.parameters.selected_onsets = pd.DataFrame({"Viewing": [None for _ in df_index]}, index=df_index)
        for i, _ in self.parameters.selected_onsets.iterrows():
            for v in self.parameters.viewings:
                try:
                    self.df_options.loc[i+(v,)]
                except KeyError:
                    continue
                self.parameters.selected_onsets.loc[i, "Viewing"] = v
                break

        if self.parameters.view_dfs:
            return self.parameters.selected_onsets

    def channel_energy_range(self, sensor, particle, particle_prefix, channel) -> tuple:
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
                self.df_grouped[sensor][particle][self.sensor_viewings(sensor)[0]][
                    particle_prefix
                ].columns
            ):
                low_energy, high_energy = self.channel_energy_range(sensor, particle, particle_prefix, channel)
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
        """Prepares the Sun-observer distance of the VDA fit (e.g. loads the SPICE kernels of the observer)"""
        self.observer.initialize_position()

    @staticmethod
    def format_timedelta(td) -> str:
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
            print(f"    Release Time : {res['Release Time']} ± {self.format_timedelta(res['Release Time Error'])}")
            print(f"    Extra Time   : {self.format_timedelta(res['Extra Time'])}")
            print(f"    APL          : {res['APL']:.2f} ± {res['APL Error']:.2f}\n")

    def compute_vda(self):
        """Fits the VDA line of each event and stores the results in self.results"""
        self.vda_fits = {}
        for index_event in self.df_options.index.unique(level=0):
            vda_points = []
            # onset times are fitted in seconds from the event start
            t0 = self.df_times.loc[index_event][self.START_TIME_COLNAME].to_pydatetime()
            t_sun_to_observer = (
                self.observer.distance(self.df_times.loc[index_event][self.START_TIME_COLNAME])
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
            self.vda_fits[index_event] = {
                "t0": t0,
                "inv_betas": inv_betas,
                "onset_seconds": onset_seconds,
                "a": a,
                "b": b,
                "b_error": b_error,
            }
