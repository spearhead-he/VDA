"""Jupyter widgets of the VDA tool.

VDA_nb_displayer creates the widgets that set the parameters of a VDA object, the background window
selection and the interactive onset selection. The plots of the widgets are created with vda_views.py.
"""
import html
from io import BytesIO

import pandas as pd
from IPython.display import Image, clear_output, display
from ipywidgets import widgets
from matplotlib import pyplot as plt

import vda_views
from vda_tool_configuration import OnsetSelection, channel_group_label


def show_figure(fig) -> None:
    """Displays a figure of vda_views in the notebook, as a .png image, and closes it"""
    buffer = BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    display(Image(buffer.getvalue()))


class VDA_nb_displayer:

    def __init__(self, vda_obj):
        ############### Widgets ###############
        self.WIDGETS_LAYOUT = widgets.Layout(width="auto")
        self.WIDGETS_STYLE = {"description_width": "initial"}

        self.vda = vda_obj

    def _bind(self, widget, parameter):
        """Sets the parameter to the value of the widget when it changes"""
        widget.observe(lambda traitlet: setattr(self.vda.parameters, parameter, traitlet["new"]), names="value")

    def display_input_file(self):
        w = widgets.Text(
            value=self.vda.parameters.input_filepath,
            placeholder="Path to .csv - Leave blank to use the datetime range below",
            description="Events file:",
            continuous_update=False,
            disabled=False,
            style=self.WIDGETS_STYLE,
            layout=self.WIDGETS_LAYOUT,
        )
        self._bind(w, "input_filepath")
        return w

    def display_date_range(self):
        """Displays the inputs of the data range, depending on the events file"""
        filepath = self.vda.parameters.input_filepath
        if not filepath:
            return self.display_date_ranges()

        try:
            file_type = self.vda.times_file_type(filepath)
        except (OSError, ValueError) as e:
            return widgets.HTML(value=f"<pre>{html.escape(str(e))}</pre>")

        if file_type == "reference times":
            return self.display_reference_hours()

        df = self.vda._read_times_file(filepath)
        with_bg = df[self.vda.BG_END_TIME_COLNAME].notna().sum()
        return widgets.Label(
            value=f"Datetime ranges from the events file: {len(df)} events, {with_bg} with a background window"
        )

    def display_reference_hours(self):
        """Hours prior to and after the reference time, for files with reference times"""
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
        self._bind(wgt_tw_prior, "bg_hours_prior")
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
        self._bind(wgt_tw_after, "bg_hours_after")
        return widgets.VBox([wgt_tw_prior, wgt_tw_after])

    def display_date_ranges(self):
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
        w = widgets.Text(
            value=self.vda.parameters.load_data_filepath,
            placeholder="Path to .pkl - Leave blank to download the data",
            description="Load the data from:",
            disabled=False,
            style=self.WIDGETS_STYLE,
            layout=self.WIDGETS_LAYOUT,
        )
        self._bind(w, "load_data_filepath")
        return w

    def display_save_data_option(self):
        w = widgets.Text(
            value=self.vda.parameters.save_data_filepath,
            placeholder="Path to .pkl - Leave blank to not save the data",
            description="Save the data to:",
            disabled=False,
            style=self.WIDGETS_STYLE,
            layout=self.WIDGETS_LAYOUT,
        )
        self._bind(w, "save_data_filepath")
        return w

    def display_particle_selection(self):
        """Channel groups, viewings and resample frequency"""
        return widgets.VBox([
            self.display_channel_groups(),
            self.display_viewings(),
            self.display_resample_frequency(),
        ])

    def _channel_options(self, sensor, species):
        """Channel numbers with their energy range, if the energies table is available"""
        channels = self.vda.parameters.AVAILABLE_CHANNELS[sensor][species]
        df_energies = getattr(self.vda, "df_energies", None)
        if df_energies is None:
            return channels
        particle_prefix = self.vda.PARTICLE_COLUMN_PREFIX[species]
        options = []
        for c in channels:
            energies = df_energies.loc[(sensor, f"{particle_prefix}_{c}")]
            options.append((f"{c}: {energies['Low Energy']:.3g}-{energies['High Energy']:.3g} MeV", c))
        return options

    def display_channel_groups(self):
        """Channel groups, with buttons to add and remove them"""
        out_options = widgets.Output()
        wrapper_channels = widgets.HBox()

        def close_options():
            out_options.clear_output(wait=False)

        def remove_channel(species, label, wgt_channel):
            channel_groups = self.vda.parameters.channel_groups
            del channel_groups[species][label]
            if len(channel_groups[species]) == 0:
                del channel_groups[species]
            wrapper_channels.children = tuple(c for c in wrapper_channels.children if c is not wgt_channel)
            wgt_channel.close()

        def set_channels(species, label, channels):
            self.vda.parameters.channel_groups[species][label]["channels"] = list(channels)

        def add_channel(sensor, species, label=None):
            """Shows the channel group with the label, or adds a new one if there is no label"""
            groups = self.vda.parameters.channel_groups.setdefault(species, {})
            if label is None:
                number = 1
                while channel_group_label(sensor, species, number) in groups:
                    number += 1
                label = channel_group_label(sensor, species, number)
                groups[label] = {"sensor": sensor, "channels": []}
            value = groups[label]["channels"]

            wgt_html = widgets.HTML(value=f"<style>p{{word-wrap: break-word; margin: 0px; text-align: center;}}</style> <p>{label}</p>")

            btn_remove = widgets.Button(description="Remove Channel", tooltip=f"Remove {label}")

            options = self._channel_options(sensor, species)
            wgt_select = widgets.SelectMultiple(options=options,
                                                value=value,
                                                description="",
                                                layout={
                                                    "min_width": "max-content",
                                                    "height": f"{2.2*len(options) + 2}ch",
                                                    "max_height": "300px"
                                                },
                                                style={"description_width": "0px"})
            wgt_select.observe(lambda traitlet: set_channels(species, label, traitlet["new"]), names="value")

            wgt_channel = widgets.VBox([wgt_html, btn_remove, wgt_select], layout={"border": "solid 1px"})
            btn_remove.on_click(lambda _: remove_channel(species, label, wgt_channel))
            wrapper_channels.children += (wgt_channel,)
            close_options()
        
        @out_options.capture(clear_output=True, wait=True)
        def show_options():
            buttons = []
            for sensor, particles in self.vda.parameters.AVAILABLE_SENSORS_PARTICLES.items():
                for species in particles:
                    button = widgets.Button(description=f"{sensor.upper()}/{species}")
                    button.on_click(lambda _, sensor=sensor, species=species: add_channel(sensor, species))
                    buttons.append(button)
            display(widgets.HBox(buttons))

        for species, groups in list(self.vda.parameters.channel_groups.items()):
            for label, group in groups.items():
                add_channel(group["sensor"], species, label)

        btn_choose = widgets.Button(description="Add Channel")
        btn_choose.on_click(lambda _: show_options())
        wrapper_btns = widgets.HBox([btn_choose, out_options])
        return widgets.VBox([wrapper_btns, wrapper_channels])

    def display_viewings(self):
        list_wgt_chk_viewings = []
        for viewing in self.vda.parameters.AVAILABLE_VIEWINGS:
            w = widgets.Checkbox(value=viewing in self.vda.parameters.viewings, 
                                 description=viewing, 
                                 disabled=False, 
                                 indent=False,
                                 layout=widgets.Layout(width="auto"))
            w.observe(lambda traitlet, viewing=viewing: self._select_viewing(viewing, traitlet["new"]),
                      names="value")
            list_wgt_chk_viewings.append(w)
        return widgets.HBox([widgets.Label("Viewings: ", style={"description_width": "max-content"})] + list_wgt_chk_viewings)

    def display_resample_frequency(self):
        wgt_resample_freq = widgets.Text(value=self.vda.parameters.resample_frequency,
                                         placeholder="Valid offset aliases string (e.g. 5min, 5T, etc) - Leave blank for no resampling",
                                         description="Resample frequency:",
                                         disabled=False,
                                         style=self.WIDGETS_STYLE,
                                         layout=self.WIDGETS_LAYOUT)
        self._bind(wgt_resample_freq, "resample_frequency")
        return wgt_resample_freq

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
        self._bind(w, "onset_method")
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

            def set_value(traitlet, parameter=parameter):
                self.vda.parameters.onset_method_parameters[parameter] = traitlet["new"]

            w.observe(set_value, names="value")

            list_param_widgets.append(w)

        return widgets.VBox(list_param_widgets)

    def display_bg_defaults(self):
        try:
            step = max(1, int(pd.Timedelta(self.vda.parameters.resample_frequency).total_seconds() // 60))
        except ValueError:
            step = 1
        df_times = getattr(self.vda, "df_times", None)
        if df_times is None:
            # before the events are created: up to one day
            longest_event = 24 * 60
        else:
            durations = df_times[self.vda.END_TIME_COLNAME] - df_times[self.vda.START_TIME_COLNAME]
            longest_event = int(durations.max().total_seconds() // 60)
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
            bg_start, bg_end = vda.bg_window(event_no)
            status = (
                f"Background {bg_start:%Y-%m-%d %H:%M} to {bg_end:%Y-%m-%d %H:%M} "
                f"({vda.bg_window_source(event_no)}), {vda.bg_window_points(event_no, times)} points"
            )
            warnings = vda.bg_window_warnings(event_no, times)
            wgt_status.value = "<br>".join([html.escape(status)] + [f"<b>{html.escape(w)}</b>" for w in warnings])
            with out_plot:
                clear_output(wait=True)
                show_figure(vda_views.plot_bg(vda, event_no))

        def show_event(_=None):
            event_no = wgt_event.value
            times = vda.df_grouped.loc[event_no].index
            bg_start, bg_end = vda.bg_window(event_no)
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
                                      ("Interactive", OnsetSelection.INTERACTIVE)], 
                             value=self.vda.parameters.onset_selection, 
                             description="Onset selection method:", 
                             disabled=False, 
                             style=self.WIDGETS_STYLE)
        self._bind(w, "onset_selection")
        return w

    def display_onset_selection(self):
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
            low_energy, high_energy = vda.channel_energy_range(sensor, particle, particle_prefix, channel)
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
                show_figure(vda_views.plot_onsets(vda, *key, selected_viewing=selected_onsets.loc[key, "Viewing"]))

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
        self._bind(w, "view_dfs")
        return w
