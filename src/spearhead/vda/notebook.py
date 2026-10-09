"""Notebook interface of the VDA tool: the parameters form, its widgets and the steps of the analysis.

The analysis itself is in analysis.py, its parameters in conf.py and the plots in views.py.
"""
import html
import os
from io import BytesIO
from warnings import filterwarnings, simplefilter

import pandas as pd
from IPython.display import Image, clear_output, display
from ipywidgets import widgets
from matplotlib import pyplot as plt
from pandas.errors import PerformanceWarning
from sunpy import log as sunpy_log

from . import views
from .analysis import VDA
from .conf import ONSET_METHOD_LABELS, OnsetSelection, VDA_parameters, channel_group_label
from .observers import OBSERVERS


def _show_figure(fig) -> None:
    """Displays a figure of views.py in the notebook, as a .png image, and closes it"""
    buffer = BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    display(Image(buffer.getvalue()))


class VDA_notebook:

    def __init__(self, vda=None):
        """Notebook interface of a VDA analysis.

        By default a new analysis with the default parameters is created. An existing VDA object can be given
        instead, e.g. to use the background or the interactive onset selection in another notebook.
        """
        simplefilter(action="ignore", category=PerformanceWarning)
        filterwarnings(action="ignore", message="Discarding nonzero nanoseconds in conversion")

        self.vda = VDA(VDA_parameters()) if vda is None else vda
        self.font_size = 22
        self.legend_font_size = 18
        self.save_plots = True
        self._apply_fonts()

        self._layout = widgets.Layout(width="auto")
        self._style = {"description_width": "initial"}

    @property
    def parameters(self):
        return self.vda.parameters

    def _apply_fonts(self):
        plt.rc("font", size=self.font_size)
        plt.rc("legend", fontsize=self.legend_font_size)

    def _show(self, df):
        if self.parameters.view_dfs:
            display(df)

    ############### Parameters form ###############
    def parameters_form(self):
        """All the parameters of the analysis, in tabs"""
        tabs = widgets.Tab(children=[
            self._events_tab(),
            self._data_tab(),
            self._channels_tab(),
            self._onsets_tab(),
            self._plots_tab(),
        ])
        for i, title in enumerate(["Events", "Data", "Energy channels", "Onsets", "Views/Plots"]):
            tabs.set_title(i, title)
        return tabs

    @staticmethod
    def _section(title, *children):
        return widgets.VBox([widgets.HTML(f"<b>{title}</b>"), *children])

    @staticmethod
    def _note(text):
        return widgets.HTML(f"<i>{html.escape(text)}</i>")

    def _bind(self, widget, parameter):
        """Sets the parameter to the value of the widget when it changes"""
        widget.observe(lambda traitlet: setattr(self.parameters, parameter, traitlet["new"]), names="value")

    def _set_enabled(self, widget, enabled):
        """Grays out (or activates again) the widget and all the widgets in it"""
        if hasattr(widget, "disabled"):
            # a single event cannot be removed
            single_event = getattr(widget, "description", "") == "Remove Event" and len(self.parameters.date_ranges) == 1
            widget.disabled = not enabled or single_event
        for child in getattr(widget, "children", ()):
            self._set_enabled(child, enabled)

    def _text_widget(self, parameter, description, placeholder):
        w = widgets.Text(
            value=getattr(self.parameters, parameter),
            placeholder=placeholder,
            description=description,
            style=self._style,
            layout=self._layout,
        )
        self._bind(w, parameter)
        return w

    ############### Events tab ###############
    def _events_tab(self):
        wgt_file = self._text_widget("input_filepath", "Events file:",
                                     "Path to .csv - Leave blank to use the datetime ranges below")
        wgt_ranges = self._date_ranges_widget()
        wrapper_file = widgets.VBox()
        shown = {"info": None}

        def update(_=None):
            filepath = self.parameters.input_filepath
            # the date ranges are grayed out as soon as an events file is typed
            self._set_enabled(wgt_ranges, not filepath)
            # the widgets of the file are only replaced when they change, not at every typed letter
            info = self._events_file_info(filepath)
            if info != shown["info"]:
                shown["info"] = info
                wrapper_file.children = self._events_file_widgets(info)

        wgt_file.observe(update, names="value")
        update()
        return widgets.VBox([
            self._section("Events file", wgt_file, wrapper_file),
            self._section("Datetime ranges (used when there is no events file)", wgt_ranges),
        ])

    def _events_file_info(self, filepath):
        """Description of the events file: ("type", text), used to show its inputs or summary"""
        if not filepath:
            return ("none", "")
        if not os.path.isfile(filepath):
            return ("note", "File not found")
        try:
            if self.vda.times_file_type(filepath) == "reference times":
                return ("reference times", "")
            df = self.vda._read_times_file(filepath)
        except ValueError as e:
            return ("error", str(e))
        with_bg = df[self.vda.BG_END_TIME_COLNAME].notna().sum()
        return ("note", f"Datetime ranges from the events file: {len(df)} events, {with_bg} with a background window")

    def _events_file_widgets(self, info):
        kind, text = info
        if kind == "none":
            return []
        if kind == "reference times":
            return [self._reference_hours_widget()]
        if kind == "error":
            return [widgets.HTML(f"<pre>{html.escape(text)}</pre>")]
        return [self._note(text)]

    def _reference_hours_widget(self):
        """Hours prior to and after the reference time, for files with reference times"""
        sliders = []
        for parameter, description in [("bg_hours_prior", "Hours prior to the reference time:"),
                                       ("bg_hours_after", "Hours after the reference time:")]:
            w = widgets.IntSlider(
                value=getattr(self.parameters, parameter),
                min=0,
                max=12,
                step=1,
                description=description,
                style=self._style,
                layout=self._layout,
            )
            self._bind(w, parameter)
            sliders.append(w)
        return widgets.VBox(sliders)

    def _date_ranges_widget(self):
        """One row with the datetime range of each event, and a button to add events"""
        date_ranges = self.parameters.date_ranges
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
                    style=self._style,
                    layout=self._layout,
                )
                wgt_start.observe(lambda traitlet, i=i: set_range_value(i, 0, traitlet["new"]), names="value")
                wgt_end = widgets.widget_datetime.NaiveDatetimePicker(
                    value=end,
                    description="end:",
                    style=self._style,
                    layout=self._layout,
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

    ############### Data tab ###############
    def _data_tab(self):
        self._wrapper_viewings = widgets.VBox([self._viewings_widget()])
        return widgets.VBox([
            self._section("Observer", self._observer_widget()),
            self._section("Viewings", self._wrapper_viewings),
            self._section("Resampling", self._text_widget(
                "resample_frequency", "Resample frequency:",
                "Valid offset aliases string (e.g. 5min, 5T, etc) - Leave blank for no resampling")),
            self._section("Advanced: load previously saved data, or save the data",
                          self._text_widget("load_data_filepath", "Load the data from:",
                                            "Path to .pkl - Leave blank to download the data"),
                          self._text_widget("save_data_filepath", "Save the data to:",
                                            "Path to .pkl - Leave blank to not save the data")),
        ])

    def _observer_widget(self):
        w = widgets.Dropdown(options=[(observer.label, name) for name, observer in OBSERVERS.items()],
                             value=self.parameters.observer,
                             description="Spacecraft:",
                             style=self._style)

        def on_change(traitlet):
            # the viewings and grouped channels are reset to the defaults of the new observer
            self.parameters.observer = traitlet["new"]
            self._wrapper_viewings.children = [self._viewings_widget()]
            self._wrapper_channels.children = [self._channels_section()]

        w.observe(on_change, names="value")
        return w

    def _viewings_widget(self):
        checkboxes = []
        for viewing in self.parameters.AVAILABLE_VIEWINGS:
            w = widgets.Checkbox(value=viewing in self.parameters.viewings,
                                 description=viewing,
                                 indent=False,
                                 layout=widgets.Layout(width="auto"))
            w.observe(lambda traitlet, viewing=viewing: self._select_viewing(viewing, traitlet["new"]),
                      names="value")
            checkboxes.append(w)
        wgt_checkboxes = widgets.HBox([widgets.Label("Viewings: ", style={"description_width": "max-content"})]
                                      + checkboxes)
        notes = []
        descriptions = self.vda.observer.VIEWING_DESCRIPTIONS
        if descriptions:
            notes.append(self._note("; ".join(f"{viewing}: {text}" for viewing, text in descriptions.items())))
        sensor_viewings = self.vda.observer.SENSOR_VIEWINGS
        if any(set(v) != set(self.parameters.AVAILABLE_VIEWINGS) for v in sensor_viewings.values()):
            text = "; ".join(f"{sensor.upper()}: {', '.join(v)}" for sensor, v in sensor_viewings.items())
            notes.append(self._note(f"Viewings of each sensor: {text}"))
        return widgets.VBox([wgt_checkboxes, *notes]) if notes else wgt_checkboxes

    def _select_viewing(self, viewing, selected):
        """Adds or removes a viewing. The checkboxes keep the order of AVAILABLE_VIEWINGS"""
        self.parameters.viewings = [
            v for v in self.parameters.AVAILABLE_VIEWINGS
            if (v == viewing and selected) or (v != viewing and v in self.parameters.viewings)
        ]

    ############### Energy channels tab ###############
    def _channels_tab(self):
        self._wrapper_channels = widgets.VBox([self._channels_section()])
        return self._wrapper_channels

    def _channels_section(self):
        # the energy ranges of the channels are the same for all the events
        start, end = self.parameters.date_ranges[0]
        self.vda.construct_energies_df(start, end)
        notes = [
            self._note(f"{self.vda.observer.label} {sensor.upper()}: the energy ranges of the channels changed on "
                       f"{date:%Y-%m-%d} ({reason}). They are read from the first event, so the events with "
                       f"{sensor.upper()} channels must all be before or all after this date.")
            for sensor, changes in self.vda.observer.ENERGY_CHANGES.items()
            for date, reason in changes
        ]
        return self._section(
            "Grouped energy channels (select multiple channels with Ctrl+click)",
            *notes,
            self._channel_groups_widget(),
        )

    def _channel_options(self, sensor, species):
        """Channel numbers with their energy range, if the energies table is available"""
        channels = self.parameters.AVAILABLE_CHANNELS[sensor][species]
        df_energies = getattr(self.vda, "df_energies", None)
        if df_energies is None:
            return channels
        particle_prefix = self.vda.PARTICLE_COLUMN_PREFIX[species]
        options = []
        for c in channels:
            energies = df_energies.loc[(sensor, f"{particle_prefix}_{c}")]
            options.append((f"{c}: {energies['Low Energy']:.3g}-{energies['High Energy']:.3g} MeV", c))
        return options

    def _channel_groups_widget(self):
        """Channel groups, with buttons to add and remove them"""
        out_options = widgets.Output()
        wrapper_channels = widgets.HBox()

        def close_options():
            out_options.clear_output(wait=False)

        def remove_channel(species, label, wgt_channel):
            channel_groups = self.parameters.channel_groups
            del channel_groups[species][label]
            if len(channel_groups[species]) == 0:
                del channel_groups[species]
            wrapper_channels.children = tuple(c for c in wrapper_channels.children if c is not wgt_channel)
            wgt_channel.close()

        def set_channels(species, label, channels):
            self.parameters.channel_groups[species][label]["channels"] = list(channels)

        def add_channel(sensor, species, label=None):
            """Shows the channel group with the label, or adds a new one if there is no label"""
            groups = self.parameters.channel_groups.setdefault(species, {})
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
            for sensor, particles in self.parameters.AVAILABLE_SENSORS_PARTICLES.items():
                for species in particles:
                    button = widgets.Button(description=f"{sensor.upper()}/{species}")
                    button.on_click(lambda _, sensor=sensor, species=species: add_channel(sensor, species))
                    buttons.append(button)
            display(widgets.HBox(buttons))

        for species, groups in list(self.parameters.channel_groups.items()):
            for label, group in groups.items():
                add_channel(group["sensor"], species, label)

        btn_choose = widgets.Button(description="Add Channel")
        btn_choose.on_click(lambda _: show_options())
        wrapper_btns = widgets.HBox([btn_choose, out_options])
        return widgets.VBox([wrapper_btns, wrapper_channels])

    ############### Onsets tab ###############
    def _onsets_tab(self):
        return widgets.VBox([
            self._section("Onset determination method", self._onset_method_widget()),
            self._section("Default background window",
                          self._bg_defaults_widget(),
                          self._note("The background window of each event can be checked and adjusted later, "
                                     "after the data are downloaded (Background windows step).")),
            self._section("Onset selection", self._onset_selection_method_widget()),
        ])

    def _onset_method_widget(self):
        """The onset method, and the parameters of the selected method"""
        w = widgets.Dropdown(options=[(ONSET_METHOD_LABELS.get(m, m), m) for m in self.parameters.AVAILABLE_ONSET_METHODS],
                             value=self.parameters.onset_method,
                             description="Onset determination method:",
                             style=self._style)
        wrapper_parameters = widgets.VBox([self._onset_method_parameters_widget()])

        def set_method(traitlet):
            # the parameters are set to the defaults of the new method
            self.parameters.onset_method = traitlet["new"]
            wrapper_parameters.children = [self._onset_method_parameters_widget()]

        w.observe(set_method, names="value")
        return widgets.VBox([w, wrapper_parameters])

    def _onset_method_parameters_widget(self):
        list_param_widgets = []
        for parameter, pinfo in self.parameters.AVAILABLE_ONSET_METHODS[self.parameters.onset_method].items():
            widget_params = {
                "value": self.parameters.onset_method_parameters.get(parameter, pinfo["default"]),
                "description": f'{parameter} | {pinfo["description"]}',
                "style": self._style,
                "layout": self._layout,
            }
            if pinfo["type"] == int:
                w = widgets.IntSlider(min=pinfo["min"], max=pinfo["max"], step=1, **widget_params)
            elif pinfo["type"] == str:
                w = widgets.Text(placeholder=pinfo["placeholder"], **widget_params)

            def set_value(traitlet, parameter=parameter):
                self.parameters.onset_method_parameters[parameter] = traitlet["new"]

            w.observe(set_value, names="value")
            list_param_widgets.append(w)

        return widgets.VBox(list_param_widgets)

    def _bg_defaults_widget(self):
        try:
            step = max(1, int(pd.Timedelta(self.parameters.resample_frequency).total_seconds() // 60))
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
            value=self.parameters.bg_after_start,
            min=0,
            max=max(longest_event, max(self.parameters.bg_after_start)),
            step=step,
            description="Default background (minutes after the start time):",
            style=self._style,
            layout=self._layout,
        )
        w.observe(
            lambda traitlet: self.vda.set_default_bg_window(*traitlet["new"]),
            names="value",
        )
        return w

    def _onset_selection_method_widget(self):
        w = widgets.Dropdown(options=[("Use all", OnsetSelection.USE_ALL),
                                      ("Interactive", OnsetSelection.INTERACTIVE)],
                             value=self.parameters.onset_selection,
                             description="Onset selection method:",
                             style=self._style)
        self._bind(w, "onset_selection")
        return w

    ############### Views/Plots tab ###############
    def _plots_tab(self):
        wgt_font = widgets.BoundedIntText(value=self.font_size, min=6, max=40, description="Font size:",
                                          style=self._style)
        wgt_legend_font = widgets.BoundedIntText(value=self.legend_font_size, min=6, max=40, description="Legend font size:",
                                                 style=self._style)
        wgt_save = widgets.Checkbox(value=self.save_plots, description="Save the VDA plots as .png files", indent=False)
        wgt_view_dfs = widgets.Checkbox(value=self.parameters.view_dfs,
                                        description="Display the produced DataFrames",
                                        indent=True,
                                        style=self._style)
        self._bind(wgt_view_dfs, "view_dfs")

        def set_font(traitlet):
            self.font_size = traitlet["new"]
            self._apply_fonts()

        def set_legend_font(traitlet):
            self.legend_font_size = traitlet["new"]
            self._apply_fonts()

        def set_save(traitlet):
            self.save_plots = traitlet["new"]

        wgt_font.observe(set_font, names="value")
        wgt_legend_font.observe(set_legend_font, names="value")
        wgt_save.observe(set_save, names="value")
        return widgets.VBox([
            self._section("Plots", wgt_font, wgt_legend_font, wgt_save),
            self._section("Tables", wgt_view_dfs),
        ])

    ############### Analysis steps ###############
    def run_data(self):
        """Creates the events, downloads (or loads) their data and groups the energy channels"""
        self.vda.construct_times_df()
        self._show(self.vda.df_times)
        self.vda.construct_particles_df()
        self.vda.group_energy_channels()
        self._show(self.vda.df_grouped)

    def background_selection(self):
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
            style=self._style,
        )
        wgt_bg = widgets.SelectionRangeSlider(
            options=[("", None)],
            description="Background:",
            continuous_update=False,
            style=self._style,
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
                _show_figure(views.plot_bg(vda, event_no))

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

    def run_onsets(self):
        """Determines the onsets and selects the ones used for the VDA (interactively, if selected)"""
        self.vda.calculate_onsets()
        self.vda.clean_onsets()
        self.vda.construct_options_df()
        self._show(self.vda.df_onsets_existing)
        self.vda.select_onsets()
        if self.parameters.onset_selection == OnsetSelection.INTERACTIVE:
            return self.onset_selection()

    def onset_selection(self):
        """Plot of the detected onsets of one channel with the selection of its viewing.

        The channels with onsets are selected with the event and channel dropdowns or the previous / next buttons.
        The selection starts from the one of vda.select_onsets().
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
            style=self._style,
        )
        wgt_channel = widgets.Dropdown(description="Channel:", style=self._style)
        wgt_viewing = widgets.Dropdown(description="Selected viewing:", style=self._style)
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
                _show_figure(views.plot_onsets(vda, *key, selected_viewing=selected_onsets.loc[key, "Viewing"]))

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

    def run_vda(self):
        """Fits the VDA line of each event, prints the results and plots them"""
        self.vda.construct_energy_channels_characteristics()
        # loading the SPICE kernels logs every frame it installs
        sunpy_level = sunpy_log.level
        sunpy_log.setLevel("WARNING")
        try:
            self.vda.define_spacecraft_parameters()
        finally:
            sunpy_log.setLevel(sunpy_level)
        self.vda.compute_vda()
        for event_no in self.vda.vda_fits:
            self.vda.print_results(event_no)
            filename = views.vda_plot_filename(self.vda, event_no) if self.save_plots else None
            _show_figure(views.plot_vda(self.vda, event_no, filename))
        self._show(self.vda.results)
