import pandas as pd
from IPython.display import display
from ipywidgets import widgets
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

    def _delete_parameter_index(self, parameter, index, cascade=False, index_sep=None):
        if index_sep is None:
            index = [index]
        else:
            index = index.split(index_sep)
        
        while True:
            par = par = self.vda.parameters.__getattribute__(parameter)
            for i in index[:-1]:
                try:
                    par = par[str(i)]
                except TypeError:
                    par = par[int(i)]
            try:
                del par[str(index[-1])]
            except TypeError:
                del par[int(index[-1])]
            
            if cascade and len(index) > 1 and len(par) == 0:
                index = index[:-1]
            else:
                break

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
            displaybox = widgets.HBox([wgt_dt_start, wgt_dt_end])
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
            displaybox = widgets.VBox([wgt_ref_times_filepath, wgt_tw_prior, wgt_tw_after])

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
            vbox = widgets.VBox([wgt_load_data_filepath])
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
            vbox = widgets.VBox([wgt_save_data, wgt_save_data_filepath])

        return vbox

    def construct_energies_df(self):
        self.vda.df_energies = pd.DataFrame({})
        df_sensors = []
        for sensor, particles in self.vda.parameters.AVAILABLE_SENSORS_PARTICLES.items():
            
            if len(particles) == 0:
                continue

            df_protons, df_electrons, energies = self.vda._epd_load(
                sensor=sensor,
                level="l2",
                startdate=self.vda.df_times.iloc[0][self.vda.START_TIME_COLNAME],
                enddate=self.vda.df_times.iloc[0][self.vda.END_TIME_COLNAME],
                viewing="sun",
                path=self.vda.DATA_PATH,
                autodownload=True,
            )
            if sensor == "het":
                flux_cols_name = "H_Flux"
                energy_bins_cols_name = "H_Bins"
            elif sensor == "ept":
                flux_cols_name = "Ion_Flux"
                energy_bins_cols_name = "Ion_Bins"
            df_protons = df_protons.rename(
                lambda x: x.replace(flux_cols_name, self.vda.PROTON_COLUMN_PREFIX),
                axis="columns",
            )
            df_electrons = df_electrons.rename(
                lambda x: x.replace("Electron_Flux", self.vda.ELECTRON_COLUMN_PREFIX),
                axis="columns",
            )

            df_energies_protons = pd.DataFrame(
                {
                    "Low Energy": energies[f"{energy_bins_cols_name}_Low_Energy"],
                    "Bin Width": energies[f"{energy_bins_cols_name}_Width"],
                },
                index=df_protons[self.vda.PROTON_COLUMN_PREFIX].columns,
            )
            df_energies_protons["High Energy"] = (
                df_energies_protons["Low Energy"] + df_energies_protons["Bin Width"]
            )

            df_energies_electrons = pd.DataFrame(
                {
                    "Low Energy": energies["Electron_Bins_Low_Energy"],
                    "Bin Width": energies["Electron_Bins_Width"],
                },
                index=df_electrons[self.vda.ELECTRON_COLUMN_PREFIX].columns,
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

        self.vda.df_energies = pd.concat(
            df_sensors,
            keys=[s for s, p in self.vda.parameters.AVAILABLE_SENSORS_PARTICLES.items() if len(p) > 0],
            names=["sensor", "channel"],
        )

        with pd.option_context("display.max_rows", None):
            display(self.vda.df_energies)

    def display_particle_selection(self):
        out_options = widgets.Output()
        wrapper_channels = widgets.HBox()
        num_channels = {}
        av_channels = self.vda.parameters.AVAILABLE_CHANNELS

        def close_options():
            out_options.clear_output(wait=False)

        def remove_channel(btn):
            key = btn.name
            self._delete_parameter_index("channel_groups", key, cascade=True, index_sep="|")
            for i, element in enumerate(wrapper_channels.children):
                btn_remove = element.children[1]
                if btn_remove.name == key:
                    element.close()
                    removed_index = i
                    break
            temp = list(wrapper_channels.children)
            del temp[removed_index]
            wrapper_channels.children = tuple(temp)

        def add_channel(btn, selected=None):
            if type(btn) == str:
                channel = btn
            else:
                channel = btn.description
            try:
                num_channels[channel] += 1
            except KeyError:
                num_channels[channel] = 1
            
            sensor, species = tuple([x.strip().lower() for x in channel.split("/")])
            label = f"{channel} Channel {num_channels[channel]}"

            value = [] if selected is None else selected
            try:
                self._change_parameter_index("channel_groups", f"{species}|{label}", {"sensor": sensor, "channels": value}, "|")
            except KeyError:
                self._change_parameter_index("channel_groups", species, {})
                self._change_parameter_index("channel_groups", f"{species}|{label}", {"sensor": sensor, "channels": value}, "|")
            
            wgt_html = widgets.HTML(value=f"<style>p{{word-wrap: break-word; margin: 0px; text-align: center;}}</style> <p>{label}</p>")

            btn_remove = widgets.Button(description="Remove Channel", tooltip=f"Remove {label}")
            btn_remove.name = f"{species}|{label}"
            btn_remove.on_click(remove_channel)

            options = av_channels[sensor][species]
            wgt_select = widgets.SelectMultiple(options=options,
                                                value=value,
                                                description=f"{species}|{label}|channels",
                                                layout={
                                                    "min_width": "max-content",
                                                    "height": f"{2.2*len(options) + 2}ch",
                                                    "max_height": "300px"
                                                },
                                                style={"description_width": "0px"})
            wgt_select.observe(lambda traitlet: self._change_parameter_index("channel_groups", 
                                                                             traitlet["owner"].description, 
                                                                             list(traitlet["new"]), "|"), 
                               names="value")

            wrapper_channels.children += (widgets.VBox([wgt_html, btn_remove, wgt_select], layout={"border": "solid 1px"}),)
            close_options()
        
        @out_options.capture(clear_output=True, wait=True)
        def show_options():
            labels = ["HET/protons", "HET/electrons", "EPT/protons", "EPT/electrons"]
            buttons = [widgets.Button(description=label) for label in labels]
            for button in buttons:
                button.on_click(add_channel)
            display(widgets.HBox(buttons))

        for species in self.vda.parameters.default_channel_groups.keys():
            for sensor, selections in self.vda.parameters.default_channel_groups[species].items():
                for selection in selections:
                    add_channel(f"{sensor}/{species}", selection)

        btn_choose = widgets.Button(description="Add Channel")
        btn_choose.on_click(lambda _: show_options())
        wrapper_btns = widgets.HBox([btn_choose, out_options])
        
        list_wgt_chk_viewings = []
        for i, viewing in enumerate(self.vda.parameters.AVAILABLE_VIEWINGS):
            w = widgets.Checkbox(value=self.vda.parameters.viewings_tt[i], 
                                 description=viewing, 
                                 disabled=False, 
                                 indent=True)
            w.observe(lambda traitlet: self._change_parameter_index("viewings_tt", 
                                                                    self.vda.parameters.AVAILABLE_VIEWINGS.index(traitlet["owner"].description), 
                                                                    traitlet["new"]),
                      names="value")
            list_wgt_chk_viewings.append(w)
        grp_viewings = widgets.HBox([widgets.Label("Viewings: ", style={"description_width": "max-content"})] + list_wgt_chk_viewings)
        

        wgt_resample_freq = widgets.Text(value=self.vda.parameters.resample_frequency,
                                         placeholder="Valid offset aliases string (e.g. 5min, 5T, etc) - Leave blank for no resampling",
                                         description="Resample frequency:",
                                         disabled=False,
                                         style=self.WIDGETS_STYLE,
                                         layout=self.WIDGETS_LAYOUT)
        wgt_resample_freq.observe(lambda traitlet: self._change_parameter("resample_frequency", 
                                                                          traitlet["new"]),
                                  names="value")
        
        display(widgets.VBox([wrapper_btns, grp_viewings, wgt_resample_freq, wrapper_channels]))

    def display_onset_method_selection(self):
        w = widgets.Dropdown(options=list(self.vda.parameters.AVAILABLE_ONSET_METHODS.keys()), 
                             value=self.vda.parameters.onset_method, 
                             description="Onset determination method:", 
                             disabled=False, 
                             style=self.WIDGETS_STYLE)
        w.observe(lambda traitlet: self._change_parameter("onset_method", 
                                                          traitlet["new"]),
                  names="value")
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

        return widgets.VBox(list_param_widgets)

        dict_wgt_onset_params = {}
        for parameter, wgt_info in AVAILABLE_ONSET_METHODS[
            wgt_onset_method.value
        ].items():
            dict_wgt_onset_params[parameter] = wgt_info["widget"](
                **wgt_info["widget_params"]
            )
        widgets.VBox(list(dict_wgt_onset_params.values()))

    def display_onset_selection_selection(self):
        w = widgets.Dropdown(options=[("Use all", 0), ("Interactive", 1), ("Custom List", 2)], 
                             value=self.vda.parameters.onset_selection, 
                             description="Onset selection method:", 
                             disabled=False, 
                             style=self.WIDGETS_STYLE)
        w.observe(lambda traitlet: self._change_parameter("onset_selection", traitlet["new"]),
                  names="value")
        return w

    def display_view_toggle(self):
        w = widgets.Checkbox(value=self.vda.parameters.view_dfs, 
                             description="Display the produced DataFrames", 
                             disabled=False, 
                             indent=True, 
                             style=self.WIDGETS_STYLE)
        w.observe(lambda traitlet: self._change_parameter("view_dfs", traitlet["new"]),
                  names="value")
        return w

    def select_onsets(self):
        temp_df = self.vda.df_options.droplevel(level=5)
        df_index = temp_df.index[~temp_df.index.duplicated(keep="first")]
        self.vda.parameters.selected_onsets = pd.DataFrame({"Viewing": [None for _ in df_index]}, index=df_index)
        if self.vda.parameters.onset_selection == 0:
            # Use all (priority by default viewings)
            for i, _ in self.vda.parameters.selected_onsets.iterrows():
                for v in self.vda.parameters.viewings:
                    try:
                        self.vda.df_options.loc[i+(v,)]
                    except KeyError:
                        continue
                    self.vda.parameters.selected_onsets.loc[i, "Viewing"] = v
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
                            display(widgets.HBox([twgt, wrb]))

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
            #                         f"Event {index_event}, {sensor}/{particle}, {self.vda.df_energies.loc[sensor, channel_low]['Low Energy']:.2f}-{self.vda.df_energies.loc[sensor, channel_high]['High Energy']:.2f} MeV, {viewing}",
            #                     )

            # assert False
            pass

    