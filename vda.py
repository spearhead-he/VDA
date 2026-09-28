import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const

from math import sqrt
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
    BG_START_TIME_COLNAME = "BG Start"
    BG_END_TIME_COLNAME = "BG End"
    END_TIME_COLNAME = "End Time"

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

    def construct_times_df(self):
        if self.parameters.input_type == 0:
            self.df_times = pd.DataFrame(
                {
                    self.BG_START_TIME_COLNAME: [self.parameters.date_start],
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
                    self.BG_START_TIME_COLNAME,
                    self.BG_END_TIME_COLNAME,
                    self.END_TIME_COLNAME
                ],
                index_col=0,
            )
            self.df_times[self.BG_START_TIME_COLNAME] = pd.to_datetime(
                self.df_times[self.BG_START_TIME_COLNAME]
            )
            self.df_times[self.BG_END_TIME_COLNAME] = pd.to_datetime(
                self.df_times[self.BG_END_TIME_COLNAME]
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
            self.df_times[self.BG_START_TIME_COLNAME] = self.df_times[
                self.REF_TIME_COLNAME
            ].apply(lambda x: x - timedelta(hours=self.parameters.bg_hours_prior))
            self.df_times[self.END_TIME_COLNAME] = self.df_times[
                self.REF_TIME_COLNAME
            ].apply(lambda x: x + timedelta(hours=self.parameters.bg_hours_after))
            self.df_times = self.df_times.drop(self.REF_TIME_COLNAME, axis="columns")

        if self.parameters.view_dfs:
            return self.df_times

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
                        startdate=row[self.BG_START_TIME_COLNAME],
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
                            (df_particle.index >= row[self.BG_START_TIME_COLNAME])
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
            for sensor, particle, viewing, particle_prefix in self._iter_sensor_particle_viewings():
                df_inner = df_event[sensor][particle][viewing][particle_prefix]
                for column_name in df_inner.columns:
                    new_kwargs = dict(kwargs)
                    new_kwargs["sensor"] = sensor
                    new_kwargs["particle"] = particle
                    new_kwargs["viewing"] = viewing
                    new_kwargs["channel"] = column_name
                    if "bg_start" in kwargs and isinstance(kwargs["bg_start"], pd.Series):
                        new_kwargs["bg_start"] = kwargs["bg_start"].loc[index_event].to_pydatetime()
                        new_kwargs["bg_end"] = kwargs["bg_end"].loc[index_event].to_pydatetime()
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

    def construct_energy_channels_characteristics(self):
        rows = []
        index = []
        for sensor, particle, particle_prefix in self._iter_sensor_particles():
            for channel in list(
                self.df_grouped[sensor][particle][self.parameters.viewings[0]][
                    particle_prefix
                ].columns
            ):
                low_energy_key = f"{particle_prefix}_{self.parameters.channel_groups[particle][channel]['channels'][0]}"
                high_energy_key = f"{particle_prefix}_{self.parameters.channel_groups[particle][channel]['channels'][-1]}"
                low_energy = self.df_energies.loc[sensor, low_energy_key]["Low Energy"]
                high_energy = self.df_energies.loc[sensor, high_energy_key]["High Energy"]
                    
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
            time_start = self.df_times.loc[index_event][self.BG_START_TIME_COLNAME].strftime("%Y-%m-%d %H:%M")
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
            t0 = self.df_times.loc[index_event][self.BG_START_TIME_COLNAME].to_pydatetime()
            t_sun_to_observer = (
                spice.get_body(
                    "Solar Orbiter",
                    self.df_times.loc[index_event][self.BG_START_TIME_COLNAME],
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

            fig, ax = plt.subplots(figsize=(10, 8))
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
            plt.title(f"Event {index_event} ({self.df_grouped.loc[index_event].index[1].to_pydatetime().strftime('%Y-%m-%d')})")
            plt.xlabel("Inverse Beta")
            plt.ylabel("Time")
            time_formatter = mdates.DateFormatter("%H:%M")
            ax.yaxis.set_major_formatter(time_formatter)
            plt.plot(
                [],
                [],
                alpha=0,
                label=f"Extra Time = {self._format_timedelta(res['Extra Time'])}",
            )
            plt.plot(
                [],
                [],
                alpha=0,
                label=f"Release Time = {res['Release Time']} +/- {self._format_timedelta(res['Release Time Error'])}",
            )
            plt.plot(
                [],
                [],
                alpha=0,
                label=f"APL = {res['APL']:.2f} +/- {res['APL Error']:.2f}",
            )
            plt.legend(bbox_to_anchor=(1, 0.6), loc="upper left")
            plt.tight_layout()
            if savefig:
                # date_str = self.df_grouped.loc[index_event].index[1].to_pydatetime().strftime('%Y-%m-%d')
                time_start_str = self.df_times.loc[index_event][self.BG_START_TIME_COLNAME].strftime("%Y-%m-%d_%H%M")
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

    def plot_bg_selection(self):
        for event_no, event in self.df_grouped.groupby(level=0):
            _, ax = plt.subplots(figsize=(10, 8))
            temp_df = event.droplevel(0)
            plt.plot(temp_df)
            bot_lim, top_lim = ax.get_ylim()
            if bot_lim <= 0:
                bot_lim = np.nanmin(temp_df.replace(0, np.nan).values)
            bg_start = self.df_times.loc[event_no][self.BG_START_TIME_COLNAME] \
                       if self.parameters.input_type == 1 \
                       else self.df_grouped.loc[event_no].index[self.parameters.onset_method_parameters["bg_start"]]
            bg_end = self.df_times.loc[event_no][self.BG_END_TIME_COLNAME] \
                     if self.parameters.input_type == 1 \
                     else self.df_grouped.loc[event_no].index[self.parameters.onset_method_parameters["bg_end"]]
            plt.fill_betweenx([0, top_lim*10],
                              bg_start,
                              bg_end,
                              color="green",
                              alpha=0.3)

            ax.set_ylabel("Flux")
            ax.set_yscale("log")
            ax.set_ylim(bot_lim/10, top_lim*10)
            ax.set_xlabel("Time")
            plt.show()
