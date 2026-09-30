import html
import pandas as pd
from IPython.display import clear_output, display
from ipywidgets import widgets
from matplotlib import pyplot as plt

from vda_tool_configuration import OnsetSelection


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

    def display_input_file(self):
        w = widgets.Text(
            value=self.vda.parameters.input_filepath,
            placeholder="Path to .csv - Leave blank to use the datetime range below",
            description="Events file:",
            disabled=False,
            style=self.WIDGETS_STYLE,
            layout=self.WIDGETS_LAYOUT,
        )
        w.observe(
            lambda traitlet: self._change_parameter("input_filepath", traitlet["new"]),
            names="value",
        )
        return w

    def display_date_range(self):
        """Displays the inputs of the data range, depending on the events file"""
        filepath = self.vda.parameters.input_filepath
        if not filepath:
            return self._display_date_ranges()

        try:
            file_type = self.vda.times_file_type(filepath)
        except (OSError, ValueError) as e:
            return widgets.HTML(value=f"<pre>{html.escape(str(e))}</pre>")

        if file_type == "reference times":
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
            return widgets.VBox([wgt_tw_prior, wgt_tw_after])

        df = self.vda._read_times_file(filepath)
        with_bg = df[self.vda.BG_END_TIME_COLNAME].notna().sum()
        return widgets.Label(
            value=f"Datetime ranges from the events file: {len(df)} events, {with_bg} with a background window"
        )

    def _display_date_ranges(self):
        """One row with the datetime range of each event, and a button to add events"""
        date_ranges = self.vda.parameters.date_ranges
        wrapper_rows = widgets.VBox()

        def set_range_value(i, j, value):
            date_range = list(date_ranges[i])
            date_range[j] = value
            date_ranges[i] = tuple(date_range)

        def remove_event(i):
            del date_ranges[i]
            show_rows()

        def add_event(_):
            date_ranges.append(date_ranges[-1])
            show_rows()

        def show_rows():
            rows = []
            for i, (start, end) in enumerate(date_ranges):
                wgt_start = widgets.widget_datetime.NaiveDatetimePicker(
                    value=start,
                    description=f"Event {i + 1} start:",
                    disabled=False,
                    style=self.WIDGETS_STYLE,
                    layout=self.WIDGETS_LAYOUT,
                )
                wgt_start.observe(lambda traitlet, i=i: set_range_value(i, 0, traitlet["new"]), names="value")
                wgt_end = widgets.widget_datetime.NaiveDatetimePicker(
                    value=end,
                    description="end:",
                    disabled=False,
                    style=self.WIDGETS_STYLE,
                    layout=self.WIDGETS_LAYOUT,
                )
                wgt_end.observe(lambda traitlet, i=i: set_range_value(i, 1, traitlet["new"]), names="value")
                btn_remove = widgets.Button(
                    description="Remove Event",
                    tooltip=f"Remove event {i + 1}",
                    disabled=len(date_ranges) == 1,
                )
                btn_remove.on_click(lambda _, i=i: remove_event(i))
                rows.append(widgets.HBox([wgt_start, wgt_end, btn_remove]))
            wrapper_rows.children = rows

        btn_add = widgets.Button(description="Add Event")
        btn_add.on_click(add_event)
        show_rows()
        return widgets.VBox([wrapper_rows, btn_add])

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
            flux_cols_name = self.vda.RAW_FLUX_COLUMN[sensor]
            energy_bins_cols_name = self.vda.RAW_ENERGY_BINS_COLUMN[sensor]
            df_protons = df_protons.rename(
                lambda x: x.replace(flux_cols_name["protons"], self.vda.PROTON_COLUMN_PREFIX),
                axis="columns",
            )
            df_electrons = df_electrons.rename(
                lambda x: x.replace(flux_cols_name["electrons"], self.vda.ELECTRON_COLUMN_PREFIX),
                axis="columns",
            )

            df_energies_protons = pd.DataFrame(
                {
                    "Low Energy": energies[f"{energy_bins_cols_name['protons']}_Low_Energy"],
                    "Bin Width": energies[f"{energy_bins_cols_name['protons']}_Width"],
                },
                index=df_protons[self.vda.PROTON_COLUMN_PREFIX].columns,
            )
            df_energies_protons["High Energy"] = (
                df_energies_protons["Low Energy"] + df_energies_protons["Bin Width"]
            )

            df_energies_electrons = pd.DataFrame(
                {
                    "Low Energy": energies[f"{energy_bins_cols_name['electrons']}_Low_Energy"],
                    "Bin Width": energies[f"{energy_bins_cols_name['electrons']}_Width"],
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
        for viewing in self.vda.parameters.AVAILABLE_VIEWINGS:
            w = widgets.Checkbox(value=viewing in self.vda.parameters.viewings, 
                                 description=viewing, 
                                 disabled=False, 
                                 indent=True)
            w.observe(lambda traitlet, viewing=viewing: self._select_viewing(viewing, traitlet["new"]),
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

    def _select_viewing(self, viewing, selected):
        """Adds or removes a viewing. The checkboxes keep the order of AVAILABLE_VIEWINGS"""
        self.vda.parameters.viewings = [
            v for v in self.vda.parameters.AVAILABLE_VIEWINGS
            if (v == viewing and selected) or (v != viewing and v in self.vda.parameters.viewings)
        ]

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

    def display_bg_defaults(self):
        try:
            step = max(1, int(pd.Timedelta(self.vda.parameters.resample_frequency).total_seconds() // 60))
        except ValueError:
            step = 1
        df_times = self.vda.df_times
        longest_event = int((df_times[self.vda.END_TIME_COLNAME] - df_times[self.vda.START_TIME_COLNAME]).max().total_seconds() // 60)
        w = widgets.IntRangeSlider(
            value=self.vda.parameters.bg_after_start,
            min=0,
            max=max(longest_event, max(self.vda.parameters.bg_after_start)),
            step=step,
            description="Default background (minutes after the start time):",
            disabled=False,
            style=self.WIDGETS_STYLE,
            layout=self.WIDGETS_LAYOUT,
        )
        w.observe(
            lambda traitlet: self.vda.set_default_bg_window(*traitlet["new"]),
            names="value",
        )
        return w

    def display_bg_selection(self):
        """Plot of the selected event with a slider to choose its background window"""
        vda = self.vda
        out_plot = widgets.Output()
        wgt_status = widgets.HTML()
        wgt_event = widgets.Dropdown(
            options=[
                (f"Event {event_no} ({vda.df_times.loc[event_no, vda.START_TIME_COLNAME]:%Y-%m-%d %H:%M})", event_no)
                for event_no in vda.df_grouped.index.unique(level=0)
            ],
            description="Event:",
            style=self.WIDGETS_STYLE,
        )
        wgt_bg = widgets.SelectionRangeSlider(
            options=[("", None)],
            description="Background:",
            continuous_update=False,
            style=self.WIDGETS_STYLE,
            layout=widgets.Layout(width="95%"),
        )
        btn_reset = widgets.Button(description="Reset to default", tooltip="Use the default background window for this event")
        # the slider is not applied to the event while it is set to the event's window
        updating = {"slider": False}

        def redraw():
            event_no = wgt_event.value
            times = vda.df_grouped.loc[event_no].index
            bg_start, bg_end = vda._bg_window(event_no)
            status = (
                f"Background {bg_start:%Y-%m-%d %H:%M} to {bg_end:%Y-%m-%d %H:%M} "
                f"({vda._bg_window_source(event_no)}), {vda._bg_window_points(event_no, times)} points"
            )
            warnings = vda._bg_window_warnings(event_no, times)
            wgt_status.value = "<br>".join([html.escape(status)] + [f"<b>{html.escape(w)}</b>" for w in warnings])
            with out_plot:
                clear_output(wait=True)
                fig, ax = plt.subplots(figsize=(12, 6))
                vda._plot_event_bg(ax, event_no)
                fig.tight_layout()
                plt.show()
                plt.close(fig)

        def show_event(_=None):
            event_no = wgt_event.value
            times = vda.df_grouped.loc[event_no].index
            bg_start, bg_end = vda._bg_window(event_no)
            # the data points closest to the event's window, inside it
            i_start = min(int(times.searchsorted(bg_start)), len(times) - 1)
            i_end = max(int(times.searchsorted(bg_end, side="right")) - 1, i_start)
            updating["slider"] = True
            wgt_bg.options = [(f"{t:%Y-%m-%d %H:%M}", t) for t in times]
            wgt_bg.index = (i_start, i_end)
            updating["slider"] = False
            redraw()

        def set_bg(traitlet):
            if updating["slider"]:
                return
            bg_start, bg_end = traitlet["new"]
            if bg_start >= bg_end:
                wgt_status.value = "<b>The background start must be before its end</b>"
                return
            vda.set_bg_window(wgt_event.value, bg_start, bg_end)
            redraw()

        def reset_bg(_):
            vda.reset_bg_window(wgt_event.value)
            show_event()

        wgt_event.observe(show_event, names="value")
        wgt_bg.observe(set_bg, names="value")
        btn_reset.on_click(reset_bg)
        show_event()
        return widgets.VBox([widgets.HBox([wgt_event, btn_reset]), wgt_bg, wgt_status, out_plot])

    def display_onset_selection_selection(self):
        w = widgets.Dropdown(options=[("Use all", OnsetSelection.USE_ALL),
                                      ("Interactive", OnsetSelection.INTERACTIVE),
                                      ("Custom List", OnsetSelection.CUSTOM_LIST)], 
                             value=self.vda.parameters.onset_selection, 
                             description="Onset selection method:", 
                             disabled=False, 
                             style=self.WIDGETS_STYLE)
        w.observe(lambda traitlet: self._change_parameter("onset_selection", traitlet["new"]),
                  names="value")
        return w

    def _display_onset_selection(self):
        """Plot of the detected onsets of one channel with the selection of its viewing.

        The channels with onsets are selected with the event and channel dropdowns or the previous / next buttons.
        """
        vda = self.vda
        selected_onsets = vda.parameters.selected_onsets
        # (event, sensor, particle, prefix, channel) of the channels with onsets
        keys = list(selected_onsets.index)
        out_plot = widgets.Output()
        wgt_status = widgets.HTML()

        def channel_label(key):
            _, sensor, particle, particle_prefix, channel = key
            low_energy, high_energy = vda._channel_energy_range(sensor, particle, particle_prefix, channel)
            return f"{sensor.upper()} {particle} {low_energy:.2f}-{high_energy:.2f} MeV"

        events = list(dict.fromkeys(key[0] for key in keys))
        wgt_event = widgets.Dropdown(
            options=[(f"Event {event_no} ({vda.df_times.loc[event_no, vda.START_TIME_COLNAME]:%Y-%m-%d %H:%M})", event_no) for event_no in events],
            description="Event:",
            style=self.WIDGETS_STYLE,
        )
        wgt_channel = widgets.Dropdown(description="Channel:", style=self.WIDGETS_STYLE)
        wgt_viewing = widgets.Dropdown(description="Selected viewing:", style=self.WIDGETS_STYLE)
        btn_previous = widgets.Button(description="◀ Previous")
        btn_next = widgets.Button(description="Next ▶")
        # changes made by the widget itself are not applied as user selections
        updating = {"flag": False}

        def current_key():
            return wgt_channel.value

        def update_status():
            key = current_key()
            event_keys = [k for k in keys if k[0] == key[0]]
            n_selected = sum(selected_onsets.loc[k, "Viewing"] is not None for k in event_keys)
            wgt_status.value = (
                f"Event {key[0]}: {n_selected} of {len(event_keys)} channels with a selected viewing"
                f" | channel {keys.index(key) + 1} of {len(keys)}"
            )
            btn_previous.disabled = keys.index(key) == 0
            btn_next.disabled = keys.index(key) == len(keys) - 1

        def redraw():
            key = current_key()
            with out_plot:
                clear_output(wait=True)
                fig = vda._plot_channel_onsets(*key, selected_viewing=selected_onsets.loc[key, "Viewing"])
                plt.show()
                plt.close(fig)

        def show_key(key):
            updating["flag"] = True
            wgt_event.value = key[0]
            wgt_channel.options = [(channel_label(k), k) for k in keys if k[0] == key[0]]
            wgt_channel.value = key
            viewings_with_onsets = [v for v in vda.parameters.viewings if v in vda.df_options.loc[key].index]
            wgt_viewing.options = [("None", None)] + [(v, v) for v in viewings_with_onsets]
            wgt_viewing.value = selected_onsets.loc[key, "Viewing"]
            updating["flag"] = False
            update_status()
            redraw()

        def on_event(traitlet):
            if not updating["flag"]:
                show_key(next(k for k in keys if k[0] == traitlet["new"]))

        def on_channel(traitlet):
            if not updating["flag"] and traitlet["new"] is not None:
                show_key(traitlet["new"])

        def on_viewing(traitlet):
            if updating["flag"]:
                return
            selected_onsets.loc[current_key(), "Viewing"] = traitlet["new"]
            update_status()
            redraw()

        def move(step):
            i = keys.index(current_key()) + step
            if 0 <= i < len(keys):
                show_key(keys[i])

        wgt_event.observe(on_event, names="value")
        wgt_channel.observe(on_channel, names="value")
        wgt_viewing.observe(on_viewing, names="value")
        btn_previous.on_click(lambda _: move(-1))
        btn_next.on_click(lambda _: move(1))
        show_key(keys[0])
        return widgets.VBox([
            widgets.HBox([wgt_event, wgt_channel, btn_previous, btn_next]),
            widgets.HBox([wgt_viewing, wgt_status]),
            out_plot,
        ])

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
        if self.vda.parameters.onset_selection in (OnsetSelection.USE_ALL, OnsetSelection.INTERACTIVE):
            # Use all (the first viewing with an onset, in the order of the viewings).
            # Interactive selection starts from the same viewings
            for i, _ in self.vda.parameters.selected_onsets.iterrows():
                for v in self.vda.parameters.viewings:
                    try:
                        self.vda.df_options.loc[i+(v,)]
                    except KeyError:
                        continue
                    self.vda.parameters.selected_onsets.loc[i, "Viewing"] = v
                    break
        if self.vda.parameters.onset_selection == OnsetSelection.INTERACTIVE:
            return self._display_onset_selection()
        elif self.vda.parameters.onset_selection == OnsetSelection.CUSTOM_LIST:
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

    