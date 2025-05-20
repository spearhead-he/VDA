import pandas as pd
from ipywidgets import widgets, VBox, HBox
from math import ceil
from matplotlib import pyplot as plt
from matplotlib import dates as mdates

from vda_tool_configuration import *


class VDA_nb_displayer:

    def __init__(self, vda_obj):
        ############### Widgets ###############
        self.WIDGETS_LAYOUT = widgets.Layout(width="auto")
        self.WIDGETS_STYLE = {"description_width": "initial"}

        self.vda = vda_obj

    def _change_parameter(self, parameter, new_value):
        self.vda.parameters.__setattr__(parameter, new_value)

    def _change_parameter_index(self, parameter, index, new_value, index_sep=None):
        if index_sep is None:
            index = [index]
        else:
            index = index.split(index_sep)
        par = self.vda.parameters.__getattribute__(parameter)
        for i in index[:-1]:
            try:
                par = par[str(i)]
            except TypeError:
                par = par[int(i)]
        try:
            par[str(index[-1])] = new_value
        except TypeError:
            par[int(index[-1])] = new_value

    def _change_parameter_df_index(self, parameter, index, col, new_value, index_sep=None):
        if index_sep is None:
            index = [index]
        else:
            index = index.split(index_sep)
        curated_indices = []
        par = self.vda.parameters.__getattribute__(parameter)
        for i in index:
            try:
                par = par.loc[str(i)]
                curated_indices.append(str(i))
            except KeyError:
                par = par.loc[int(i)]
                curated_indices.append(int(i))
        self.vda.parameters.__getattribute__(parameter).loc[tuple(curated_indices), col] = new_value

    def display_input_type(self):
        w = widgets.Dropdown(
            options=[
                ("Custom datetime range (1 event)", 0),
                ("File with datetime ranges", 1),
                ("File with reference datetimes", 2),
            ],
            value=self.vda.parameters.input_type,
            description="How the events will be provided:",
            disabled=False,
            style=self.WIDGETS_STYLE,
        )
        w.observe(
            lambda traitlet: self._change_parameter("input_type", traitlet["new"]),
            names="value",
        )
        return w

    def display_date_range(self):
        if self.vda.parameters.input_type == 0:
            wgt_dt_start = widgets.widget_datetime.NaiveDatetimePicker(
                value=self.vda.parameters.date_start,
                description="Datetime range start:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_dt_start.observe(
                lambda traitlet: self._change_parameter("date_start", traitlet["new"]),
                names="value",
            )
            wgt_dt_end = widgets.widget_datetime.NaiveDatetimePicker(
                value=self.vda.parameters.date_end,
                description="Datetime range end:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_dt_end.observe(
                lambda traitlet: self._change_parameter("date_end", traitlet["new"]),
                names="value",
            )
            displaybox = HBox([wgt_dt_start, wgt_dt_end])
        elif self.vda.parameters.input_type == 1:
            wgt_dt_range_filepath = widgets.Text(
                value=self.vda.parameters.date_range_filepath,
                placeholder="Path to .csv",
                description="Datetime range file:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_dt_range_filepath.observe(
                lambda traitlet: self._change_parameter(
                    "date_range_filepath", traitlet["new"]
                ),
                names="value",
            )
            displaybox = wgt_dt_range_filepath
        elif self.vda.parameters.input_type == 2:
            wgt_ref_times_filepath = widgets.Text(
                value=self.vda.parameters.reference_times_filepath,
                placeholder="Path to .csv",
                description="Reference times file:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_ref_times_filepath.observe(
                lambda traitlet: self._change_parameter(
                    "reference_times_filepath", traitlet["new"]
                ),
                names="value",
            )
            wgt_tw_prior = widgets.IntSlider(
                value=self.vda.parameters.bg_hours_prior,
                min=0,
                max=12,
                step=1,
                description="Hours prior to the reference time:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_tw_prior.observe(
                lambda traitlet: self._change_parameter(
                    "bg_hours_prior", traitlet["new"]
                ),
                names="value",
            )
            wgt_tw_after = widgets.IntSlider(
                value=self.vda.parameters.bg_hours_after,
                min=0,
                max=12,
                step=1,
                description="Hours after the reference time:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_tw_after.observe(
                lambda traitlet: self._change_parameter(
                    "bg_hours_after", traitlet["new"]
                ),
                names="value",
            )
            displaybox = VBox([wgt_ref_times_filepath, wgt_tw_prior, wgt_tw_after])

        return displaybox

    def display_load_data_option(self):
        w = widgets.Checkbox(
            value=self.vda.parameters.load_data,
            description="Load data",
            disabled=False,
            indent=True,
        )
        w.observe(
            lambda traitlet: self._change_parameter("load_data", traitlet["new"]),
            names="value",
        )
        return w

    def display_save_data_option(self):
        if self.vda.parameters.load_data:
            wgt_load_data_filepath = widgets.Text(
                value=self.vda.parameters.load_data_filepath,
                placeholder="Path to .pkl",
                description="File with saved DataFrame:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_load_data_filepath.observe(
                lambda traitlet: self._change_parameter(
                    "load_data_filepath", traitlet["new"]
                ),
                names="value",
            )
            vbox = VBox([wgt_load_data_filepath])
        else:
            wgt_save_data = widgets.Checkbox(
                value=self.vda.parameters.save_data,
                description="Save downloaded data",
                disabled=False,
                indent=True,
            )
            wgt_save_data.observe(
                lambda traitlet: self._change_parameter("save_data", traitlet["new"]),
                names="value",
            )
            wgt_save_data_filepath = widgets.Text(
                value=self.vda.parameters.save_data_filepath,
                placeholder="Path with .pkl extension",
                description="File to save data DataFrame:",
                disabled=False,
                style=self.WIDGETS_STYLE,
                layout=self.WIDGETS_LAYOUT,
            )
            wgt_save_data_filepath.observe(
                lambda traitlet: self._change_parameter(
                    "save_data_filepath", traitlet["new"]
                ),
                names="value",
            )
            vbox = VBox([wgt_save_data, wgt_save_data_filepath])

        return vbox

    def display_particle_selection(self):
        list_wgt_chk_sensors_particles = []
        for sensor, particles in self.vda.parameters.AVAILABLE_SENSORS_PARTICLES.items():
            for i, particle in enumerate(particles):
                w = widgets.Checkbox(
                    value=self.vda.parameters.sensors_particles_tt[sensor][i],
                    description=f"{sensor}/{particle}",
                    disabled=False,
                    indent=True,
                )
                w.observe(
                    lambda traitlet: self._change_parameter_index(
                        "sensors_particles_tt",
                        f'{(parts := traitlet["owner"].description.split("/"))[0]}/{self.vda.parameters.AVAILABLE_SENSORS_PARTICLES[parts[0]].index(parts[1])}',
                        traitlet["new"],
                        "/",
                    ),
                    names="value",
                )
                list_wgt_chk_sensors_particles.append(w)
        grp_sensors_particles = VBox(list_wgt_chk_sensors_particles)

        # for i, sensor in enumerate(self.vda.parameters.AVAILABLE_SENSORS):
        #     w = widgets.Checkbox(
        #         value=self.vda.parameters.sensors_tt[i],
        #         description=sensor,
        #         disabled=False,
        #         indent=True,
        #     )
        #     w.observe(
        #         lambda traitlet: self._change_parameter_index(
        #             "sensors_tt",
        #             self.vda.parameters.AVAILABLE_SENSORS.index(
        #                 traitlet["owner"].description
        #             ),
        #             traitlet["new"],
        #         ),
        #         names="value",
        #     )
        #     list_wgt_chk_sensors.append(w)
        # grp_sensors = VBox(list_wgt_chk_sensors)

        # list_wgt_chk_particles = []
        # for i, particle in enumerate(self.vda.parameters.AVAILABLE_PARTICLES):
        #     w = widgets.Checkbox(
        #         value=self.vda.parameters.particles_tt[i],
        #         description=particle,
        #         disabled=False,
        #         indent=True,
        #     )
        #     w.observe(
        #         lambda traitlet: self._change_parameter_index(
        #             "particles_tt",
        #             self.vda.parameters.AVAILABLE_PARTICLES.index(
        #                 traitlet["owner"].description
        #             ),
        #             traitlet["new"],
        #         ),
        #         names="value",
        #     )
        #     list_wgt_chk_particles.append(w)
        # grp_particles = VBox(list_wgt_chk_particles)

        list_wgt_chk_viewings = []
        for i, viewing in enumerate(self.vda.parameters.AVAILABLE_VIEWINGS):
            w = widgets.Checkbox(
                value=self.vda.parameters.viewings_tt[i],
                description=viewing,
                disabled=False,
                indent=True,
            )
            w.observe(
                lambda traitlet: self._change_parameter_index(
                    "viewings_tt",
                    self.vda.parameters.AVAILABLE_VIEWINGS.index(
                        traitlet["owner"].description
                    ),
                    traitlet["new"],
                ),
                names="value",
            )
            list_wgt_chk_viewings.append(w)
        grp_viewings = VBox(list_wgt_chk_viewings)

        wgt_resample_freq = widgets.Text(
            value=self.vda.parameters.resample_frequency,
            placeholder="Valid offset aliases string (e.g. 5min, 5T, etc) - Leave blank for no resampling",
            description="Resample frequency:",
            disabled=False,
            style=self.WIDGETS_STYLE,
            layout=self.WIDGETS_LAYOUT,
        )
        wgt_resample_freq.observe(
            lambda traitlet: self._change_parameter(
                "resample_frequency", traitlet["new"]
            ),
            names="value",
        )

        return VBox([HBox([grp_sensors_particles, grp_viewings]), wgt_resample_freq])

    def display_groupings(self):
        list_wgt_group_size = []
        for sensor, particles in self.vda.parameters.sensors_particles.items():
            for particle in particles:
                w = widgets.IntSlider(
                    value=self.vda.parameters.group_sizes[sensor][particle],
                    min=1,
                    max=6,
                    step=1,
                    description=f"{sensor}/{particle}",
                    disabled=False,
                    style=self.WIDGETS_STYLE,
                    layout=self.WIDGETS_LAYOUT,
                )
                w.observe(
                    lambda traitlet: self._change_parameter_index(
                        "group_sizes",
                        traitlet["owner"].description,
                        traitlet["new"],
                        "/"
                    ),
                    names="value",
                )
                list_wgt_group_size.append(w)

        return VBox(list_wgt_group_size)

    def display_onset_method_selection(self):
        w = widgets.Dropdown(
            options=list(self.vda.parameters.AVAILABLE_ONSET_METHODS.keys()),
            value=self.vda.parameters.onset_method,
            description="Onset determination method:",
            disabled=False,
            style=self.WIDGETS_STYLE,
        )
        w.observe(
            lambda traitlet: self._change_parameter("onset_method", traitlet["new"]),
            names="value",
        )
        return w

    def display_onset_method_parameters(self):
        list_param_widgets = []
        for parameter, pinfo in self.vda.parameters.AVAILABLE_ONSET_METHODS[
            self.vda.parameters.onset_method
        ].items():
            widget_type = None
            widget_params = {}
            widget_params["disabled"] = False
            widget_params["style"] = self.WIDGETS_STYLE
            widget_params["layout"] = self.WIDGETS_LAYOUT
            if pinfo["type"] == int:
                widget_type = widgets.IntSlider
                widget_params["value"] = pinfo["default"]
                widget_params["min"] = pinfo["min"]
                widget_params["max"] = pinfo["max"]
                widget_params["step"] = 1
                widget_params["description"] = f'{parameter} | {pinfo["description"]}'
            elif pinfo["type"] == str:
                widget_type = widgets.Text
                widget_params["value"] = pinfo["default"]
                widget_params["placeholder"] = pinfo["placeholder"]
                widget_params["description"] = f'{parameter} | {pinfo["description"]}'

            w = widget_type(**widget_params)
            w.observe(
                lambda traitlet: self._change_parameter_index(
                    "onset_method_parameters",
                    traitlet["owner"].description.split("|")[0].strip(),
                    traitlet["new"],
                ),
                names="value",
            )

            list_param_widgets.append(w)

        return VBox(list_param_widgets)

        dict_wgt_onset_params = {}
        for parameter, wgt_info in AVAILABLE_ONSET_METHODS[
            wgt_onset_method.value
        ].items():
            dict_wgt_onset_params[parameter] = wgt_info["widget"](
                **wgt_info["widget_params"]
            )
        VBox(list(dict_wgt_onset_params.values()))

    def display_onset_selection_selection(self):
        w = widgets.Dropdown(
            options=[("Use all", 0), ("Interactive", 1), ("Custom List", 2)],
            value=self.vda.parameters.onset_selection,
            description="Onset selection method:",
            disabled=False,
            style=self.WIDGETS_STYLE,
        )
        w.observe(
            lambda traitlet: self._change_parameter("onset_selection", traitlet["new"]),
            names="value",
        )
        return w

    def display_view_toggle(self):
        w = widgets.Checkbox(
            value=self.vda.parameters.view_dfs,
            description="Display the produced DataFrames",
            disabled=False,
            indent=True,
            style=self.WIDGETS_STYLE,
        )
        w.observe(
            lambda traitlet: self._change_parameter("view_dfs", traitlet["new"]),
            names="value",
        )
        return w

    def select_onsets(self):
        temp_df = self.vda.df_options.droplevel(level=5)
        df_index = temp_df.index[~temp_df.index.duplicated(keep="first")]
        self.vda.parameters.selected_onsets = pd.DataFrame({"Viewing": [None for _ in df_index]}, index=df_index)
        if self.vda.parameters.onset_selection == 0:
            # Use all (priority by default viewings)
            for i, row in self.vda.parameters.selected_onsets.iterrows():
                for v in self.vda.parameters.viewings:
                    try:
                        self.vda.df_options.loc[i+(v,)]
                    except KeyError:
                        continue
                    row["Viewing"] = v
        elif self.vda.parameters.onset_selection == 1:
            time_formatter = mdates.DateFormatter("%H:%M")
            for event_no, event in self.vda.df_grouped.groupby(level=0):
                temp_df = event.droplevel(0)
                for sensor, particles in self.vda.parameters.sensors_particles.items():
                    for particle in particles:
                        if particle == "protons":
                            particle_prefix = self.vda.PROTON_COLUMN_PREFIX
                        elif particle == "electrons":
                            particle_prefix = self.vda.ELECTRON_COLUMN_PREFIX
                        columns = temp_df[sensor][particle][self.vda.parameters.viewings[0]][particle_prefix].columns
                        for column in columns:
                            onset_found = False
                            nplots = len(self.vda.parameters.viewings)
                            ncols = 3
                            if nplots <= ncols:
                                nrows = 1
                                ncols = nplots
                            else:
                                nrows = ceil(nplots/ncols)
                            fig, axs = plt.subplots(nrows,
                                                    ncols,
                                                    figsize=(14, 8))
                            try:
                                axs_flat = axs.flatten()
                            except AttributeError:
                                axs_flat = [axs]
                            for ax in axs_flat[len(self.vda.parameters.viewings):]:
                                ax.axis("off")
                            for ax, viewing in zip(axs_flat, self.vda.parameters.viewings):
                                ax.set_title(viewing)
                                try:
                                    onset_results = self.vda.df_onsets_existing.loc[(event_no,
                                                                                sensor,
                                                                                particle,
                                                                                viewing,
                                                                                particle_prefix,
                                                                                column)]
                                    onset_found = True
                                except KeyError:
                                    continue

                                ax.plot(temp_df[sensor][particle][viewing][particle_prefix][column].fillna(0).ffill(), label="Data")
                                xlim = ax.get_xlim()
                                ylim = ax.get_ylim()
                                ax.fill_betweenx([0, ylim[1]],
                                                onset_results["Background Start"],
                                                onset_results["Background End"],
                                                color="green",
                                                alpha=0.3,
                                                label="BG Sample")
                                ax.hlines(onset_results["Method Specific"]["bg_level"],
                                        xlim[0],
                                        xlim[1],
                                        color="green",
                                        linestyles="dashed",
                                        label=f'BG ({onset_results["Method Specific"]["bg_level"]:.2f})')
                                ax.hlines(onset_results["Method Specific"]["threshold"],
                                        xlim[0],
                                        xlim[1],
                                        color="red",
                                        linestyles="dashed",
                                        label=f'Threshold ({onset_results["Method Specific"]["threshold"]:.2f})')
                                ax.vlines(onset_results["Onset Time"],
                                        0,
                                        ylim[1],
                                        color="purple",
                                        linestyles="dashed",
                                        label=f'Onset ({onset_results["Onset Time"].strftime("%H:%M")})')

                                ax.set_xlim(xlim)
                                ax.xaxis.set_major_formatter(time_formatter)
                                # ax.set_ylim(top=ylim[1])
                                ax.set_yscale("log")

                                ax.legend()
                            
                            if not onset_found:
                                plt.close()
                                continue

                            if "-" in column:
                                channels = column.split("-")
                                low_energy_key = channels[0]
                                high_energy_key = channels[1]
                            else:
                                low_energy_key = column
                                high_energy_key = column
                            energy_range_str = f"{self.vda.df_energies.loc[(sensor, low_energy_key), 'Low Energy']:.2f}-{self.vda.df_energies.loc[(sensor, high_energy_key), 'High Energy']:.2f}"
                            
                            twgt = widgets.Label(
                                value=f"Event {event_no} ({temp_df.index[0].to_pydatetime().strftime('%Y-%m-%d')}) | {sensor}/{particle} ({energy_range_str} MeV):",
                                style=self.WIDGETS_STYLE,
                                layout=self.WIDGETS_LAYOUT
                            )
                            wrb = widgets.RadioButtons(
                                options=[None] + [v for v in self.vda.parameters.viewings 
                                                    if v in self.vda.df_options.loc[(event_no,
                                                                                sensor,
                                                                                particle,
                                                                                particle_prefix,
                                                                                column)].index],
                                index=0,
                                description=f"{event_no}/{sensor}/{particle}/{particle_prefix}/{column}",
                                orientation="horizontal",
                                style=self.WIDGETS_STYLE,
                                layout=self.WIDGETS_LAYOUT
                            )
                            wrb.observe(
                                lambda traitlet: self._change_parameter_df_index(
                                    "selected_onsets",
                                    traitlet["owner"].description,
                                    "Viewing",
                                    traitlet["new"],
                                    "/",
                                ),
                                names="value",
                            )
                            display(HBox([twgt, wrb]))

                            plt.suptitle(f"Detected onsets for event {event_no} ({temp_df.index[0].to_pydatetime().strftime('%Y-%m-%d')}) | {sensor}/{particle} ({energy_range_str} MeV)")
                            plt.tight_layout()
                            plt.show()
        elif self.vda.parameters.onset_selection == 2:
            # # Custom list
            # df_selections = pd.DataFrame({})
            # for index_event, df_event in self.df_options.groupby(level=0):
            #     for sensor, particles in self.vda.parameters.sensors_particles.items():
            #         for particle in particles:
            #             if particle == "protons":
            #                 particle_prefix = self.PROTON_COLUMN_PREFIX
            #             elif particle == "electrons":
            #                 particle_prefix = self.ELECTRON_COLUMN_PREFIX
            #             for channel, df_channel in df_event.loc[
            #                 index_event, sensor, particle, particle_prefix
            #             ].groupby(level=0):
            #                 channel_low = (channels := channel.split("-"))[0]
            #                 channel_high = channels[1]
            #                 for viewing, df_viewing in df_channel.groupby(level=1):
            #                     self._plot_onset(
            #                         self.df_grouped.loc[index_event][
            #                             sensor,
            #                             particle,
            #                             viewing,
            #                             particle_prefix,
            #                             channel,
            #                         ],
            #                         df_channel.loc[channel, viewing][
            #                             "Onset Time"
            #                         ].to_pydatetime(),
            #                         df_channel.loc[channel, viewing][
            #                             "Background Start"
            #                         ].to_pydatetime(),
            #                         df_channel.loc[channel, viewing][
            #                             "Background End"
            #                         ].to_pydatetime(),
            #                         f"Event {index_event}, {sensor}/{particle}, {self.df_energies.loc[sensor, channel_low]['Low Energy']:.2f}-{self.df_energies.loc[sensor, channel_high]['High Energy']:.2f} MeV, {viewing}",
            #                     )

            # assert False
            pass

    