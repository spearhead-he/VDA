import astrospice
import numpy as np
import pandas as pd
import astropy.units as u

from math import sqrt
from os import getcwd
from datetime import timezone, datetime, timedelta

from matplotlib import pyplot as plt
from matplotlib import dates as mdates
from sunpy.coordinates import HeliocentricInertial
from solo_epd_loader import epd_load
from pyonset import Onset, BootstrapWindow


class VDA:

    def __init__(self, parameters):
        self.parameters = parameters

    ############### Reference Times DF ###############
    @property
    def EVENT_INDEX_NAME(self):
        return "Event No"

    @property
    def REF_TIME_COLNAME(self):
        return "Reference Time"

    @property
    def START_TIME_COLNAME(self):
        return "Start Time"

    @property
    def END_TIME_COLNAME(self):
        return "End Time"

    ############### Particle Data ###############
    @property
    def DATA_PATH(self):
        return f"{getcwd()}/particle_data"

    @property
    def PROTON_COLUMN_PREFIX(self):
        return "H_Flux"

    @property
    def ELECTRON_COLUMN_PREFIX(self):
        return "Electron_Flux"

    ############### Onset Selection ###############
    @property
    def VIEWINGS_HIERARCHY(self):
        return ["sun", "north", "south", "asun", "omni"]

    ############### VDA ###############
    @property
    def C(self):
        return 299_792_458

    @property
    def AU_TO_M_RATIO(self):
        return 1.495978707e11

    @property
    def M_REST(self):
        return {"protons": 938.27, "electrons": 0.511}

    def construct_times_df(self):
        if self.parameters.input_type == 0:
            self.df_times = pd.DataFrame(
                {
                    self.START_TIME_COLNAME: [self.parameters.date_start],
                    self.END_TIME_COLNAME: [self.parameters.date_end],
                },
                index=[1],
            )
        elif self.parameters.input_type == 1:
            self.df_times = pd.read_csv(
                self.parameters.date_range_filepath,
                sep=",",
                header=0,
                names=[
                    self.EVENT_INDEX_NAME,
                    self.START_TIME_COLNAME,
                    self.END_TIME_COLNAME,
                ],
                index_col=0,
            )
            self.df_times[self.START_TIME_COLNAME] = pd.to_datetime(
                self.df_times[self.START_TIME_COLNAME]
            )
            self.df_times[self.END_TIME_COLNAME] = pd.to_datetime(
                self.df_times[self.END_TIME_COLNAME]
            )
        elif self.parameters.input_type == 2:
            self.df_times = pd.read_csv(
                self.parameters.reference_times_filepath,
                sep=",",
                header=0,
                names=[self.EVENT_INDEX_NAME, self.REF_TIME_COLNAME],
                index_col=0,
            )
            self.df_times[self.REF_TIME_COLNAME] = pd.to_datetime(
                self.df_times[self.REF_TIME_COLNAME]
            )
            self.df_times[self.START_TIME_COLNAME] = self.df_times[
                self.REF_TIME_COLNAME
            ].apply(lambda x: x - timedelta(hours=self.parameters.bg_hours_prior))
            self.df_times[self.END_TIME_COLNAME] = self.df_times[
                self.REF_TIME_COLNAME
            ].apply(lambda x: x + timedelta(hours=self.parameters.bg_hours_after))
            self.df_times = self.df_times.drop(self.REF_TIME_COLNAME, axis="columns")

        # self.df_times = self.df_times.map(lambda x: x.replace(tzinfo=timezone.utc))

        if self.parameters.view_dfs:
            return self.df_times

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
                    df_protons, df_electrons, _ = epd_load(
                        sensor=sensor,
                        level="l2",
                        startdate=row[self.START_TIME_COLNAME],
                        enddate=row[self.END_TIME_COLNAME],
                        viewing=viewing,
                        path=self.DATA_PATH,
                        autodownload=True,
                    )
                    if "protons" in particles:
                        if sensor == "het":
                            flux_cols_name = "H_Flux"
                        elif sensor == "ept":
                            flux_cols_name = "Ion_Flux"
                        df_protons = df_protons[
                            [c for c in df_protons.columns if c[0] == flux_cols_name]
                        ]
                        # df_protons.index = df_protons.index.tz_localize(timezone.utc)
                        df_protons = df_protons[
                            (df_protons.index >= row[self.START_TIME_COLNAME])
                            & (df_protons.index <= row[self.END_TIME_COLNAME])
                        ]
                        if (
                            self.parameters.resample_frequency is not None
                            and self.parameters.resample_frequency != ""
                        ):
                            df_protons = df_protons.resample(
                                self.parameters.resample_frequency, origin="start"
                            ).mean()
                            df_protons.index = df_protons.index.floor("min")
                        df_protons = pd.concat(
                            [df_protons],
                            keys=[(sensor, "protons", viewing)],
                            axis="columns",
                        )
                        df_protons = df_protons.rename(
                            lambda x: x.replace(
                                flux_cols_name, self.PROTON_COLUMN_PREFIX
                            ),
                            axis="columns",
                        )
                        df_row = pd.concat([df_row, df_protons], axis="columns")
                    if "electrons" in particles:
                        if sensor == "het" or sensor == "ept":
                            flux_cols_name = "Electron_Flux"
                        df_electrons = df_electrons[
                            [c for c in df_electrons.columns if c[0] == flux_cols_name]
                        ]
                        # df_electrons.index = df_electrons.index.tz_localize(
                        #     timezone.utc
                        # )
                        df_electrons = df_electrons[
                            (df_electrons.index >= row[self.START_TIME_COLNAME])
                            & (df_electrons.index <= row[self.END_TIME_COLNAME])
                        ]
                        if (
                            self.parameters.resample_frequency is not None
                            and self.parameters.resample_frequency != ""
                        ):
                            df_electrons = df_electrons.resample(
                                self.parameters.resample_frequency, origin="start"
                            ).mean()
                            df_electrons.index = df_electrons.index.floor("min")
                        df_electrons = pd.concat(
                            [df_electrons],
                            keys=[(sensor, "electrons", viewing)],
                            axis="columns",
                        )
                        df_electrons = df_electrons.rename(
                            lambda x: x.replace(
                                flux_cols_name, self.ELECTRON_COLUMN_PREFIX
                            ),
                            axis="columns",
                        )
                        df_row = pd.concat([df_row, df_electrons], axis="columns")
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

    def construct_energies_df(self):
        self.df_energies = pd.DataFrame({})
        df_sensors = []
        for sensor, particles in self.parameters.sensors_particles.items():
            
            if len(particles) == 0:
                continue

            df_protons, df_electrons, energies = epd_load(
                sensor=sensor,
                level="l2",
                startdate=self.df_times.iloc[0][self.START_TIME_COLNAME],
                enddate=self.df_times.iloc[0][self.END_TIME_COLNAME],
                viewing="sun",
                path=self.DATA_PATH,
                autodownload=True,
            )
            if sensor == "het":
                flux_cols_name = "H_Flux"
                energy_bins_cols_name = "H_Bins"
            elif sensor == "ept":
                flux_cols_name = "Ion_Flux"
                energy_bins_cols_name = "Ion_Bins"
            df_protons = df_protons.rename(
                lambda x: x.replace(flux_cols_name, self.PROTON_COLUMN_PREFIX),
                axis="columns",
            )
            df_electrons = df_electrons.rename(
                lambda x: x.replace("Electron_Flux", self.ELECTRON_COLUMN_PREFIX),
                axis="columns",
            )

            df_energies_protons = pd.DataFrame(
                {
                    "Low Energy": energies[f"{energy_bins_cols_name}_Low_Energy"],
                    "Bin Width": energies[f"{energy_bins_cols_name}_Width"],
                },
                index=df_protons[self.PROTON_COLUMN_PREFIX].columns,
            )
            df_energies_protons["High Energy"] = (
                df_energies_protons["Low Energy"] + df_energies_protons["Bin Width"]
            )

            df_energies_electrons = pd.DataFrame(
                {
                    "Low Energy": energies["Electron_Bins_Low_Energy"],
                    "Bin Width": energies["Electron_Bins_Width"],
                },
                index=df_electrons[self.ELECTRON_COLUMN_PREFIX].columns,
            )
            df_energies_electrons["High Energy"] = (
                df_energies_electrons["Low Energy"] + df_energies_electrons["Bin Width"]
            )

            if "protons" in particles and "electrons" in particles:
                df_sensors.append(pd.concat([df_energies_protons, df_energies_electrons]))
            elif "protons" in particles:
                df_sensors.append(pd.concat([df_energies_protons]))
            elif "electrons" in particles:
                df_sensors.append(pd.concat([df_energies_electrons]))

        self.df_energies = pd.concat(
            df_sensors,
            keys=[s for s, p in self.parameters.sensors_particles.items() if len(p) > 0],
            names=["sensor", "channel"],
        )

        if self.parameters.view_dfs:
            return self.df_energies

    def _group_channels_de(
        self,
        df_to_group: pd.DataFrame,
        energy_bins_width: list[float],
        column_prefix: str = "",
        group_size: int = 2,
    ) -> pd.DataFrame:
        # I = ΣI_n*ΔE_n / ΣΔE_n
        
        if group_size == 1:
            return df_to_group
        
        channels = list(range(len(df_to_group.columns)))
        grouped_channels = [
            channels[c : c + group_size] for c in range(0, len(channels), group_size)
        ]
        grouped_all = {}
        for group in grouped_channels:
            de = sum([energy_bins_width[group[i]] for i, _ in enumerate(group)])
            grouped_series = df_to_group.iloc[:, group[0]] * energy_bins_width[group[0]]
            for i, _ in enumerate(group[1:]):
                grouped_series = grouped_series.add(
                    df_to_group.iloc[:, group[i]] * energy_bins_width[group[i]],
                    fill_value=0,
                )
            grouped_all[
                f"{column_prefix}{'_' if column_prefix != '' else ''}{group[0]}-{column_prefix}{'_' if column_prefix != '' else ''}{group[-1]}"
            ] = (grouped_series / de)
        df_grouped = pd.DataFrame(grouped_all)
        return df_grouped

    def group_energy_channels(self):
        self.df_grouped = pd.DataFrame({})
        for sensor, particles in self.parameters.sensors_particles.items():
            for particle in particles:
                if particle == "protons":
                    particle_prefix = self.PROTON_COLUMN_PREFIX
                elif particle == "electrons":
                    particle_prefix = self.ELECTRON_COLUMN_PREFIX
                for viewing in self.parameters.viewings:
                    df_temp = self._group_channels_de(
                        self.df_data[sensor][particle][viewing][particle_prefix],
                        list(self.df_energies.loc[sensor]["Bin Width"]),
                        column_prefix=particle_prefix,
                        group_size=self.parameters.group_sizes[sensor][particle],
                    )
                    df_temp = pd.concat(
                        [df_temp],
                        keys=[(sensor, particle, viewing, particle_prefix)],
                        axis="columns",
                    )
                    self.df_grouped = pd.concat(
                        [self.df_grouped, df_temp], axis="columns"
                    )

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
        if type(bg_start) is int:
            bg_start = series.index[bg_start]
        if type(bg_end) is int:
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
        if type(bg_start) is int:
            bg_start = series.index[bg_start]
        if type(bg_end) is int:
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
    ) -> dict:
        df_onsets = pd.DataFrame({})
        for index_event, df_event in df.groupby(level=0):
            for sensor, particles in self.parameters.sensors_particles.items():
                for particle in particles:
                    if particle == "protons":
                        particle_prefix = self.PROTON_COLUMN_PREFIX
                    elif particle == "electrons":
                        particle_prefix = self.ELECTRON_COLUMN_PREFIX
                    for viewing in self.parameters.viewings:
                        for column_name in (
                            df_inner := df_event[sensor][particle][viewing][
                                particle_prefix
                            ]
                        ).columns:
                            kwargs["sensor"] = sensor
                            kwargs["particle"] = particle
                            kwargs["viewing"] = viewing
                            kwargs["channel"] = column_name
                            try:
                                onset_time, bg_start, bg_stop, method_specific = (
                                    self._onset_detection(
                                        df_inner[column_name].droplevel(
                                            0, axis="index"
                                        ),
                                        method,
                                        **kwargs,
                                    )
                                )
                                df_onsets = pd.concat(
                                    [
                                        df_onsets,
                                        pd.DataFrame(
                                            {
                                                "Onset Time": [onset_time],
                                                "Background Start": [bg_start],
                                                "Background End": [bg_stop],
                                                "Method Specific": [method_specific],
                                            },
                                            index=[
                                                [index_event],
                                                [sensor],
                                                [particle],
                                                [viewing],
                                                [particle_prefix],
                                                [column_name],
                                            ],
                                        ),
                                    ]
                                )
                            except Exception as e:
                                print(index_event, type(e).__name__, kwargs)
                                df_onsets = pd.concat(
                                    [
                                        df_onsets,
                                        pd.DataFrame(
                                            {
                                                "Onset Time": [pd.NaT],
                                                "Background Start": [pd.NaT],
                                                "Background End": [pd.NaT],
                                                "Method Specific": [None],
                                            },
                                            index=[
                                                [index_event],
                                                [sensor],
                                                [particle],
                                                [viewing],
                                                [particle_prefix],
                                                [column_name],
                                            ],
                                        ),
                                    ]
                                )
        df_onsets.index.names = [
            self.EVENT_INDEX_NAME,
            "sensor",
            "particle",
            "viewing",
            "prefix",
            "channels",
        ]
        return df_onsets

    def calculate_onsets(self):
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

        # if self.parameters.view_dfs:
        #     return self.df_options

    def _plot_onset(
        self,
        series: pd.Series,
        onset_time: datetime,
        bg_start_time: datetime,
        bg_end_time: datetime,
        title: str,
        vlines: dict = None,
        hlines: dict = None,
    ) -> None:
        ax = series.fillna(0).plot(title=title, logy=True, label="Data")
        ax.set_ylim((ylim := ax.get_ylim())[0] * 0.01, ylim[1] * 100)
        ax.axvline(onset_time, linestyle="--", label="Onset time")
        ylim_top = ax.get_ylim()[1]
        ax.fill_between(
            [bg_start_time, bg_end_time],
            0,
            ylim_top,
            color="green",
            alpha=0.25,
            label="BG sample",
        )
        if vlines is not None:
            for label, line_info in vlines.items():
                ax.axvline(
                    line_info["value"], label=label, **line_info.get("lineargs", {})
                )
        if hlines is not None:
            for label, line_info in hlines.items():
                ax.axhline(
                    line_info["value"], label=label, **line_info.get("lineargs", {})
                )
        ax.set_ylim(top=ylim_top)
        ax.legend()
        plt.axes(ax)
        plt.tight_layout()
        plt.show()

    def construct_energy_channels_characteristics(self):
        self.df_channels_chars = pd.DataFrame({})
        for sensor, particles in self.parameters.sensors_particles.items():
            for particle in particles:
                if particle == "protons":
                    particle_prefix = self.PROTON_COLUMN_PREFIX
                elif particle == "electrons":
                    particle_prefix = self.ELECTRON_COLUMN_PREFIX
                for channel in list(
                    self.df_grouped[sensor][particle][self.parameters.viewings[0]][
                        particle_prefix
                    ].columns
                ):
                    low_energy_key = channel
                    high_energy_key = channel
                    if "-" in channel:
                        channels = channel.split("-")
                        low_energy_key = channels[0]
                        high_energy_key = channels[1]
                    low_energy = self.df_energies.loc[sensor, low_energy_key]["Low Energy"]
                    high_energy = self.df_energies.loc[sensor, high_energy_key][
                        "High Energy"
                    ]
                    
                    geo_mean = sqrt(low_energy) * sqrt(high_energy)
                    inv_beta = 1 / sqrt(
                        1 - (1 / (1 + geo_mean / self.M_REST[particle])) ** 2
                    )
                    self.df_channels_chars = pd.concat(
                        [
                            self.df_channels_chars,
                            pd.DataFrame(
                                {
                                    "Geomagnetic Mean": [geo_mean],
                                    "Inverse Beta": [inv_beta],
                                },
                                index=[[sensor], [particle], [channel]],
                            ),
                        ]
                    )
        self.df_channels_chars.index.names = ["sensor", "particle", "channel"]

        if self.parameters.view_dfs:
            return self.df_channels_chars

    def define_spacecraft_parameters(self):
        kernels = astrospice.registry.get_kernels("solar orbiter", "predict")
        solo_kernel = kernels[0]
        coverage = solo_kernel.coverage("SOLAR ORBITER")
        print(coverage.iso)

    def plot(self, savefig: bool = True):
        heliocentric = HeliocentricInertial()
        for index_event, df_event in self.df_options.groupby(level=0):
            vda_points = []
            t_sun_to_observer = (
                astrospice.generate_coords(
                    "SOLAR ORBITER",
                    self.df_times.loc[index_event][self.START_TIME_COLNAME],
                )
                .transform_to(heliocentric)
                .distance.to(u.au)[0]
                .value
                * self.AU_TO_M_RATIO
                / self.C
            )
            for i, row in self.parameters.selected_onsets.loc[index_event].iterrows():
                if row["Viewing"] is None:
                    continue
                sensor, particle, particle_prefix, channel = i
                vda_points.append(
                    (self.df_channels_chars.loc[
                        sensor, particle, channel
                    ]["Inverse Beta"],
                    self.df_onsets_existing.loc[
                        index_event,
                        sensor,
                        particle,
                        row["Viewing"],
                        particle_prefix,
                        channel
                    ]["Onset Time"]
                    .to_pydatetime()
                    .timestamp())
                )

            if len(vda_points) < 2:
                # Not enough points for the linear regression
                # Consider throughing warning
                print(f"Not enough onset points in event {index_event}.")
                continue

            vda_points = sorted(vda_points, key=lambda x: x[0])
            inv_betas = np.array([p[0] for p in vda_points])
            timestamps = np.array([p[1] for p in vda_points])

            try:
                p, V = np.polyfit(inv_betas, timestamps, 1, cov=True)
                a = p[0]
                b = p[1]
                a_error = np.sqrt(V[0][0])
                b_error = np.sqrt(V[1][1])
            except ValueError:
                # not enough points for cov matrix
                print(
                    f"Not enough points for covariance matrix generation in event {index_event}"
                )
                a, b = np.polyfit(inv_betas, timestamps, 1)
                a_error = 0
                b_error = 0

            fig, ax = plt.subplots(figsize=(10, 8))
            ax.scatter(
                inv_betas,
                [datetime.fromtimestamp(t) for t in timestamps],
                color="black",
            )
            ax.plot(
                inv_betas,
                [datetime.fromtimestamp(a * x + b) for x in inv_betas],
                label="Linear Regression",
                color="blue",
            )
            ax.fill_between(
                inv_betas,
                [datetime.fromtimestamp(a * x + b - 2 * b_error) for x in inv_betas],
                [datetime.fromtimestamp(a * x + b + 2 * b_error) for x in inv_betas],
                color="blue",
                alpha=0.1,
            )
            plt.title(f"Event {index_event} ({self.df_grouped.loc[index_event].index[1].to_pydatetime().strftime('%Y-%m-%d')})")
            plt.xlabel("Inverse Beta")
            plt.ylabel("Time")
            time_formatter = mdates.DateFormatter("%H:%M")
            ax.yaxis.set_major_formatter(time_formatter)
            plt.plot(
                [],
                [],
                alpha=0,
                label=f"Extra Time = {str(timedelta(seconds=t_sun_to_observer)).split('.')[0]}",
            )
            plt.plot(
                [],
                [],
                alpha=0,
                label=f"Release Time = {datetime.fromtimestamp(b + t_sun_to_observer).strftime('%Y-%m-%d %H:%M:%S')} +/- {str(timedelta(seconds=b_error)).split('.')[0]}",
            )
            plt.plot(
                [],
                [],
                alpha=0,
                label=f"APL = {a / t_sun_to_observer:.2f} +/- {a_error / t_sun_to_observer:.2f}",
            )
            plt.legend(bbox_to_anchor=(1, 0.6), loc="upper left")
            plt.tight_layout()
            if savefig:
                # date_str = self.df_grouped.loc[index_event].index[1].to_pydatetime().strftime('%Y-%m-%d')
                time_start_str = self.df_times.loc[index_event]["Start Time"].strftime("%Y-%m-%d_%H%M")
                time_end_str = self.df_times.loc[index_event]["End Time"].strftime("%Y-%m-%d_%H%M")
                date_str = f"{time_start_str}_{time_end_str}"
                particles_str = "_".join([f'{s}-{p}' for s, ps in self.parameters.sensors_particles.items() for p in ps])
                freq_str = self.parameters.resample_frequency if self.parameters.resample_frequency != "" else "noresample"
                filename = f"{date_str}_{particles_str}_{freq_str}.png"
                plt.savefig(filename)
            plt.show()

    def plot_bg_selection(self):
        for event_no, event in self.df_grouped.groupby(level=0):
            _, ax = plt.subplots(figsize=(10, 8))
            temp_df = event.droplevel(0)
            plt.plot(temp_df)
            bot_lim, top_lim = ax.get_ylim()
            if bot_lim <= 0:
                bot_lim = np.nanmin(temp_df.replace(0, np.nan).values)
            plt.fill_betweenx(
                [0, top_lim*10],
                self.df_grouped.loc[event_no].index[self.parameters.onset_method_parameters["bg_start"]],
                self.df_grouped.loc[event_no].index[self.parameters.onset_method_parameters["bg_end"]],
                color="green",
                alpha=0.3)

            ax.set_ylabel("Flux")
            ax.set_yscale("log")
            ax.set_ylim(bot_lim/10, top_lim*10)
            ax.set_xlabel("Time")
            plt.show()
