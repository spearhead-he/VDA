"""Notebook interfaces of the VDA tool: the parameters form, its widgets and the steps of the analysis.

VDA_notebook is the interface of the VDA of one or more events from one spacecraft (vda_multi_event.ipynb),
VDA_multi_spacecraft_notebook the one of the VDA of one event from several spacecraft (vda_multi_spacecraft.ipynb),
and VDA_custom_points_notebook the one of the VDA of onset points given by the user (vda_custom_points.ipynb).
The analysis itself is in analysis.py, its parameters in conf.py and the plots in views.py.
"""
import html
import os
from copy import deepcopy
from datetime import datetime
from io import BytesIO
from warnings import filterwarnings, simplefilter

import pandas as pd
from IPython.display import Image, clear_output, display
from ipywidgets import widgets
from matplotlib import pyplot as plt
from pandas.errors import PerformanceWarning
from sunpy import log as sunpy_log

from . import views
from .analysis import VDA, CustomPointsVDA, PARTICLE_NAMES, points_df, read_points_file
from .conf import DEFAULT_DATA_PATH, ONSET_METHOD_LABELS, OnsetSelection, VDA_parameters, channel_group_label
from .observers import OBSERVERS


def _show_figure(fig) -> None:
    """Displays a figure of views.py in the notebook, as a .png image, and closes it"""
    buffer = BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    display(Image(buffer.getvalue()))


def _is_connection_error(error) -> bool:
    """Whether the error is a failed connection to a data server (e.g. aborted, refused, SSL error or timeout)"""
    return any(cls.__name__ in ("ConnectionError", "SSLError", "Timeout", "TimeoutError") for cls in type(error).__mro__)


def _print_results(label, res, format_timedelta) -> None:
    """Prints the results of a VDA fit (a row of VDA.results) under the label"""
    print(label)
    if pd.isna(res.get("APL")):
        print("    No results (not enough onset points)\n")
        return
    print(f"    Release Time : {res['Release Time']} ± {format_timedelta(res['Release Time Error'])}")
    print(f"    Extra Time   : {format_timedelta(res['Extra Time'])}")
    print(f"    APL          : {res['APL']:.2f} ± {res['APL Error']:.2f}\n")


def _compute_vda(vda) -> None:
    """VDA fit of the events of the analysis"""
    vda.construct_energy_channels_characteristics()
    # loading the SPICE kernels logs every frame it installs, and JPL Horizons every location it gives
    sunpy_level = sunpy_log.level
    sunpy_log.setLevel("WARNING")
    try:
        vda.define_spacecraft_parameters()
        vda.compute_vda()
    finally:
        sunpy_log.setLevel(sunpy_level)


class _NotebookBase:
    """The parts shared by the notebook interfaces: the fonts, saving and reference times of the plots"""

    def __init__(self):
        simplefilter(action="ignore", category=PerformanceWarning)
        filterwarnings(action="ignore", message="Discarding nonzero nanoseconds in conversion")

        self.font_size = 22
        self.legend_font_size = 18
        self.save_plots = True
        # (label, datetime) reference times shown in the VDA plots, e.g. SXR peak, type II and III radio bursts
        self.reference_times = []
        self._apply_fonts()

        self._layout = widgets.Layout(width="auto")
        self._style = {"description_width": "initial"}

    def _apply_fonts(self):
        plt.rc("font", size=self.font_size)
        plt.rc("legend", fontsize=self.legend_font_size)

    @staticmethod
    def _section(title, *children):
        return widgets.VBox([widgets.HTML(f"<b>{title}</b>"), *children])

    @staticmethod
    def _note(text):
        return widgets.HTML(f"<i>{html.escape(text)}</i>")

    @staticmethod
    def _error_text(text):
        return f'<span style="color: red"><b>Error:</b> {html.escape(text)}</span>'

    ############### Views/Plots tab ###############
    # where the reference times are shown, in the title of their section
    _REFERENCE_TIMES_PLOTS = "in the VDA plots"
    _SAVE_PLOTS_DESCRIPTION = "Save the VDA plots as .png files"

    def _plots_tab(self):
        return widgets.VBox(self._plots_sections())

    def _plots_sections(self):
        """Font sizes, saving of the plots and reference times"""
        wgt_font = widgets.BoundedIntText(value=self.font_size, min=6, max=40, description="Font size:",
                                          style=self._style)
        wgt_legend_font = widgets.BoundedIntText(value=self.legend_font_size, min=6, max=40, description="Legend font size:",
                                                 style=self._style)
        wgt_save = widgets.Checkbox(value=self.save_plots, description=self._SAVE_PLOTS_DESCRIPTION, indent=False)

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
        return [
            self._section("Plots", wgt_font, wgt_legend_font, wgt_save),
            self._section(f"Reference times (shown at inverse beta 0 {self._REFERENCE_TIMES_PLOTS})",
                          self._reference_times_widget()),
        ]

    def _reference_times_widget(self):
        """One row with the label and datetime of each reference time, and a button to add one"""
        references = self.reference_times
        wrapper_rows = widgets.VBox()

        def set_reference_value(i, j, value):
            reference = list(references[i])
            reference[j] = value
            references[i] = tuple(reference)

        def remove_reference(i):
            del references[i]
            show_rows()

        def add_reference(_):
            references.append((f"Reference {len(references) + 1}", self._default_reference_time()))
            show_rows()

        def show_rows():
            rows = []
            for i, (label, time) in enumerate(references):
                wgt_label = widgets.Text(value=label, description="Label:", style=self._style)
                wgt_label.observe(lambda traitlet, i=i: set_reference_value(i, 0, traitlet["new"]), names="value")
                wgt_time = widgets.widget_datetime.NaiveDatetimePicker(
                    value=time,
                    description="time:",
                    style=self._style,
                    layout=self._layout,
                )
                wgt_time.observe(lambda traitlet, i=i: set_reference_value(i, 1, traitlet["new"]), names="value")
                btn_remove = widgets.Button(description="Remove", tooltip=f"Remove {label}")
                btn_remove.on_click(lambda _, i=i: remove_reference(i))
                rows.append(widgets.HBox([wgt_label, wgt_time, btn_remove]))
            wrapper_rows.children = rows

        btn_add = widgets.Button(description="Add Reference Time")
        btn_add.on_click(add_reference)
        show_rows()
        return widgets.VBox([wrapper_rows, btn_add])

    def _default_reference_time(self):
        """Datetime of a new reference time, to be changed"""
        return datetime.now().replace(second=0, microsecond=0)


class VDA_notebook(_NotebookBase):

    def __init__(self, vda=None):
        """Notebook interface of a VDA analysis.

        By default a new analysis with the default parameters is created. An existing VDA object can be given
        instead, e.g. to use the background or the interactive onset selection in another notebook.
        """
        self.vda = VDA(VDA_parameters()) if vda is None else vda
        super().__init__()

    @property
    def parameters(self):
        return self.vda.parameters

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
            self._section("Downloaded data", self._data_path_widget()),
            self._section("Advanced: load previously saved data, or save the data",
                          self._text_widget("load_data_filepath", "Load the data from:",
                                            "Path to .pkl - Leave blank to download the data"),
                          self._text_widget("save_data_filepath", "Save the data to:",
                                            "Path to .pkl - Leave blank to not save the data")),
        ])

    def _data_path_widget(self) -> widgets.Text:
        return self._text_widget("data_path", "Data folder:",
                                 f"Folder of the downloaded data files - Leave blank for {DEFAULT_DATA_PATH} "
                                 "in the current directory")

    def _observer_widget(self):
        w = widgets.Dropdown(options=[(observer.label, name) for name, observer in OBSERVERS.items()],
                             value=self.parameters.observer,
                             description="Spacecraft:",
                             style=self._style)

        # error of the data of the observer, e.g. no data for the first event
        self._observer_error = widgets.HTML()

        def on_change(traitlet):
            # the viewings and grouped channels are reset to the defaults of the new observer
            self.parameters.observer = traitlet["new"]
            self._wrapper_viewings.children = [self._viewings_widget()]
            self._wrapper_channels.children = [self._channels_section()]

        w.observe(on_change, names="value")
        return widgets.VBox([w, self._observer_error])

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
        # without data, the channels are listed without their energy ranges
        self.vda.df_energies = None
        error = ""
        try:
            self.vda.construct_energies_df(start, end)
        except Exception as e:
            self.vda.df_energies = None
            if _is_connection_error(e):
                advice = "the data server may be unavailable, read the channels again later."
            else:
                advice = "choose another spacecraft, or change the events and read the channels again."
            error = self._error_text(
                f"{e}. The energy ranges of the {self.vda.observer.label} channels could not be read from the data "
                f"of the first event: {advice}"
            )
        self._observer_error.value = error
        notes = [
            self._note(f"{self.vda.observer.label} {sensor.upper()}: the energy ranges of the channels changed on "
                       f"{date:%Y-%m-%d} ({reason}). They are read from the first event, so the events with "
                       f"{sensor.upper()} channels must all be before or all after this date.")
            for sensor, changes in self.vda.observer.ENERGY_CHANGES.items()
            for date, reason in changes
        ]
        wgt_channel_groups = self._channel_groups_widget()
        # without the data of the observer the grouped channels cannot be changed
        self._set_enabled(wgt_channel_groups, not error)
        wgt_error = []
        if error:
            btn_retry = widgets.Button(description="Read the channels again",
                                       tooltip="Read the energy ranges of the channels from the first event",
                                       layout=widgets.Layout(width="max-content"))
            btn_retry.on_click(lambda _: setattr(self._wrapper_channels, "children", [self._channels_section()]))
            wgt_error = [widgets.HBox([widgets.HTML(error), btn_retry])]
        return self._section(
            "Grouped energy channels (select multiple channels with Ctrl+click)",
            *wgt_error,
            *notes,
            wgt_channel_groups,
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
            key = (sensor, f"{particle_prefix}_{c}")
            if key not in df_energies.index:
                # a channel without energy range in the data
                options.append((str(c), c))
                continue
            energies = df_energies.loc[key]
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
    _REFERENCE_TIMES_PLOTS = "in the VDA plot of the event they are in"

    def _plots_tab(self):
        wgt_view_dfs = widgets.Checkbox(value=self.parameters.view_dfs,
                                        description="Display the produced DataFrames",
                                        indent=True,
                                        style=self._style)
        self._bind(wgt_view_dfs, "view_dfs")
        return widgets.VBox([*self._plots_sections(), self._section("Tables", wgt_view_dfs)])

    def _default_reference_time(self):
        # at the start of the first event
        return self.parameters.date_ranges[0][0]

    def _event_references(self, vda, event_no):
        """The reference times inside the data range of the event"""
        start, end = vda.df_times.loc[event_no, [vda.START_TIME_COLNAME, vda.END_TIME_COLNAME]]
        return [(label, time) for label, time in self.reference_times
                if time is not None and start <= pd.Timestamp(time) <= end]

    ############### Analysis steps ###############
    def run_data(self):
        """Creates the events, downloads (or loads) their data and groups the energy channels"""
        self.vda.construct_times_df()
        self._show(self.vda.df_times)
        if getattr(self.vda, "df_energies", None) is None:
            # not read by the parameters form, e.g. without data for its spacecraft and first event
            self.vda.construct_energies_df()
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

    def _determine_onsets(self):
        self.vda.calculate_onsets()
        self.vda.clean_onsets()
        self.vda.construct_options_df()
        self._show(self.vda.df_onsets_existing)
        self.vda.select_onsets()

    def run_onsets(self):
        """Determines the onsets and selects the ones used for the VDA (interactively, if selected)"""
        self._determine_onsets()
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
        if not keys:
            return widgets.HTML("No onsets were determined")
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
        _compute_vda(self.vda)
        for event_no in self.vda.vda_fits:
            self.vda.print_results(event_no)
            filename = views.vda_plot_filename(self.vda, event_no) if self.save_plots else None
            _show_figure(views.plot_vda(self.vda, event_no, filename, self._event_references(self.vda, event_no)))
        self._show(self.vda.results)


# Spacecraft of the multi-spacecraft notebook by default, with data and results for the default event
DEFAULT_SPACECRAFT = ("solo", "sta", "soho")

# Parameters of VDA_multi_spacecraft_notebook shared by all the spacecraft, in the order they are copied
SHARED_PARAMETERS = ("date_ranges", "bg_after_start", "data_path", "onset_method", "onset_method_parameters",
                     "onset_selection", "view_dfs")


class VDA_multi_spacecraft_notebook(VDA_notebook):

    def __init__(self, spacecraft=DEFAULT_SPACECRAFT):
        """Notebook interface of the VDA of one event from several spacecraft.

        Each spacecraft has its own analysis (a VDA_notebook with its VDA), with its viewings, resample frequency and
        grouped channels. The event, the onset determination and selection, and the plots are shared by all the
        spacecraft: they are the parameters of self.vda, copied to the analysis of each spacecraft when the steps run.
        """
        super().__init__()
        self.show_spacecraft_plots = False
        # the analysis of each spacecraft, kept when it is unselected
        self._notebooks = {}
        self._settings = {}
        self.spacecraft = []
        for name in spacecraft:
            self.select_spacecraft(name)
        # the analyses of the spacecraft with data, after run_data()
        self.analyses = {}
        self.results = None

    def notebook(self, name) -> VDA_notebook:
        """The analysis of the spacecraft (e.g. tool.notebook("sta").parameters for its parameters)"""
        if name not in self._notebooks:
            self._notebooks[name] = VDA_notebook(VDA(VDA_parameters(observer=name)))
            # the fonts of the plots are the ones of this notebook
            self._apply_fonts()
        return self._notebooks[name]

    def select_spacecraft(self, name, selected=True):
        """Adds or removes a spacecraft. The spacecraft keep the order of OBSERVERS"""
        if name not in OBSERVERS:
            raise ValueError(f"Unknown observer {name!r}. Use one of: {', '.join(OBSERVERS)}")
        self.notebook(name)
        self.spacecraft = [n for n in OBSERVERS if (n == name and selected) or (n != name and n in self.spacecraft)]

    def _sync(self, notebook):
        """Copies the shared parameters to the analysis of a spacecraft"""
        for parameter in SHARED_PARAMETERS:
            setattr(notebook.parameters, parameter, deepcopy(getattr(self.parameters, parameter)))

    @staticmethod
    def _label(name):
        return OBSERVERS[name].label

    ############### Parameters form ###############
    def parameters_form(self):
        """All the parameters of the analysis, in tabs"""
        tabs = widgets.Tab(children=[
            self._event_tab(),
            self._spacecraft_tab(),
            self._onsets_tab(),
            self._plots_tab(),
        ])
        for i, title in enumerate(["Event", "Spacecraft", "Onsets", "Views/Plots"]):
            tabs.set_title(i, title)
        return tabs

    def _event_tab(self):
        start, end = self.parameters.date_ranges[0]
        wgt_start = widgets.widget_datetime.NaiveDatetimePicker(value=start, description="Start:",
                                                                style=self._style, layout=self._layout)
        wgt_end = widgets.widget_datetime.NaiveDatetimePicker(value=end, description="end:",
                                                              style=self._style, layout=self._layout)

        def set_range(_):
            self.parameters.date_ranges = [(wgt_start.value, wgt_end.value)]

        wgt_start.observe(set_range, names="value")
        wgt_end.observe(set_range, names="value")
        return widgets.VBox([
            self._section("Datetime range of the event (the same for all the spacecraft)",
                          widgets.HBox([wgt_start, wgt_end])),
            self._section("Downloaded data (the same for all the spacecraft)", self._data_path_widget()),
        ])

    def _spacecraft_tab(self):
        tabs = widgets.Tab()

        def show_settings():
            tabs.children = [self._spacecraft_settings(name) for name in self.spacecraft]
            for i, name in enumerate(self.spacecraft):
                tabs.set_title(i, self._label(name))

        def select(name, selected):
            self.select_spacecraft(name, selected)
            show_settings()

        checkboxes = []
        for name in OBSERVERS:
            w = widgets.Checkbox(value=name in self.spacecraft, description=self._label(name), indent=False,
                                 layout=widgets.Layout(width="auto"))
            w.observe(lambda traitlet, name=name: select(name, traitlet["new"]), names="value")
            checkboxes.append(w)
        show_settings()
        return widgets.VBox([
            self._section("Spacecraft", widgets.HBox(checkboxes)),
            self._section("Data and energy channels of each spacecraft", tabs),
        ])

    def _spacecraft_settings(self, name):
        """Viewings, resample frequency and grouped channels of the spacecraft (created when it is first shown)"""
        if name not in self._settings:
            notebook = self.notebook(name)
            # the energy ranges of the channels are read from the data of the event
            self._sync(notebook)
            notebook._observer_error = widgets.HTML()
            notebook._wrapper_viewings = widgets.VBox([notebook._viewings_widget()])
            self._settings[name] = widgets.VBox([
                notebook._observer_error,
                self._section("Viewings", notebook._wrapper_viewings),
                self._section("Resampling", notebook._text_widget(
                    "resample_frequency", "Resample frequency:",
                    "Valid offset aliases string (e.g. 5min, 5T, etc) - Leave blank for no resampling")),
                notebook._channels_tab(),
            ])
        return self._settings[name]

    _REFERENCE_TIMES_PLOTS = "in the combined plot and the VDA plots of the spacecraft"

    def _plots_tab(self):
        tab = super()._plots_tab()
        w = widgets.Checkbox(value=self.show_spacecraft_plots,
                             description="Show the VDA plot of each spacecraft, besides the combined plot",
                             indent=False)
        w.observe(lambda traitlet: setattr(self, "show_spacecraft_plots", traitlet["new"]), names="value")
        section_plots = tab.children[0]
        section_plots.children += (w,)
        return tab

    ############### Analysis steps ###############
    def _check_analyses(self):
        if not self.analyses:
            raise ValueError("No spacecraft with data: run the data step with at least one spacecraft with data")

    def run_data(self):
        """Creates the event, downloads the data of each spacecraft and groups their energy channels.

        A spacecraft without data is left out of the next steps, and its error is printed.
        """
        if not self.spacecraft:
            raise ValueError("No spacecraft: select at least one in the Spacecraft tab")
        self.analyses = {}
        for name in self.spacecraft:
            notebook = self.notebook(name)
            self._sync(notebook)
            # the energy ranges of the channels are read again, from the data of the event
            notebook.vda.df_energies = None
            print(f"{self._label(name)}:")
            try:
                notebook.run_data()
            except Exception as e:
                if _is_connection_error(e):
                    print("    Left out, the data could not be downloaded (the data server may be unavailable, "
                          f"run the data step again later): {e}")
                else:
                    print(f"    Left out, no data: {e}")
                continue
            self.analyses[name] = notebook
        self._check_analyses()

    def _per_spacecraft(self, widget_of):
        """The widget of the spacecraft selected in a dropdown, created when it is first shown"""
        self._check_analyses()
        wgt_spacecraft = widgets.Dropdown(options=[(self._label(name), name) for name in self.analyses],
                                          description="Spacecraft:", style=self._style)
        wrapper = widgets.VBox()
        shown = {}

        def show(_=None):
            name = wgt_spacecraft.value
            if name not in shown:
                shown[name] = widget_of(self.analyses[name])
            wrapper.children = [shown[name]]

        wgt_spacecraft.observe(show, names="value")
        show()
        return widgets.VBox([wgt_spacecraft, wrapper])

    def background_selection(self):
        """Plot of the event with a slider to choose its background window, for the selected spacecraft"""
        return self._per_spacecraft(lambda notebook: notebook.background_selection())

    def run_onsets(self):
        """Determines the onsets of each spacecraft and selects the ones used for the VDA (interactively, if selected)"""
        self._check_analyses()
        for name, notebook in self.analyses.items():
            if self.parameters.view_dfs:
                print(f"{self._label(name)}:")
            notebook._determine_onsets()
        if self.parameters.onset_selection == OnsetSelection.INTERACTIVE:
            return self._per_spacecraft(lambda notebook: notebook.onset_selection())

    def run_vda(self):
        """Fits the VDA line of each spacecraft, prints the results and plots them together"""
        self._check_analyses()
        vdas = {}
        for name, notebook in self.analyses.items():
            _compute_vda(notebook.vda)
            vdas[self._label(name)] = notebook.vda
        self.results = pd.DataFrame(
            [vda.results.loc[1] if 1 in vda.results.index else pd.Series(dtype=object) for vda in vdas.values()],
            index=pd.Index(list(vdas), name="Spacecraft"),
        )
        for label, vda in vdas.items():
            _print_results(label, self.results.loc[label], vda.format_timedelta)
        references = self._event_references(next(iter(vdas.values())), 1)
        if self.show_spacecraft_plots:
            for name, notebook in self.analyses.items():
                if 1 in notebook.vda.vda_fits:
                    filename = f"{name}_{views.vda_plot_filename(notebook.vda, 1)}" if self.save_plots else None
                    _show_figure(views.plot_vda(notebook.vda, 1, filename, references))
        filename = views.vda_combined_plot_filename(vdas) if self.save_plots else None
        _show_figure(views.plot_vda_combined(vdas, 1, filename, references))
        self._show(self.results)


class VDA_custom_points_notebook(_NotebookBase):
    _REFERENCE_TIMES_PLOTS = "in the VDA plot"
    _SAVE_PLOTS_DESCRIPTION = "Save the VDA plot as a .png file"

    def __init__(self):
        """Notebook interface of the VDA fit of onset points given by the user.

        The points are in groups (e.g. one per spacecraft), each with a label, the distance of its observer from the
        Sun (AU) and its points: (onset time, energy, particle), with the energy as a single value or a range in MeV
        (e.g. "13.5" or "13-16"). They are entered in the form, loaded into it with load_points(), or read from the
        points file (points_filepath) instead of the form.
        """
        super().__init__()
        # points file used instead of the groups of the form, if set
        self.points_filepath = ""
        # groups of the form: {"label", "distance", "points": [[onset time, energy, particle], ...]}
        self.groups = [self._new_group(1)]
        self.title = ""
        # the VDA fit of each group, after run_vda()
        self.vdas = {}
        self.results = None

    @staticmethod
    def _new_group(number):
        return {"label": f"Group {number}", "distance": 1.0, "points": []}

    def load_points(self, filepath):
        """Loads the points of a points file into the groups of the form, to see or change them"""
        self.groups = [
            {
                "label": label,
                "distance": vda.distance,
                "points": [
                    [f"{row['Onset Time']:%Y-%m-%d %H:%M:%S}",
                     f"{row['Low Energy']:g}" if row["Low Energy"] == row["High Energy"]
                     else f"{row['Low Energy']:g}-{row['High Energy']:g}",
                     row["Particle"]]
                    for _, row in vda.points.iterrows()
                ],
            }
            for label, vda in read_points_file(filepath).items()
        ]

    ############### Parameters form ###############
    def parameters_form(self):
        """All the parameters, in tabs"""
        tabs = widgets.Tab(children=[self._points_tab(), self._plots_tab()])
        for i, title in enumerate(["Points", "Views/Plots"]):
            tabs.set_title(i, title)
        return tabs

    def _set_enabled(self, widget, enabled):
        """Grays out (or activates again) the widget and all the widgets in it"""
        if hasattr(widget, "disabled"):
            # a single group cannot be removed
            single_group = getattr(widget, "description", "") == "Remove Group" and len(self.groups) == 1
            widget.disabled = not enabled or single_group
        for child in getattr(widget, "children", ()):
            self._set_enabled(child, enabled)

    def _points_tab(self):
        wgt_file = widgets.Text(value=self.points_filepath, description="Points file:", style=self._style,
                                layout=self._layout,
                                placeholder="Path to .csv - Leave blank to use the groups of points below")
        wgt_file_info = widgets.HTML()
        wgt_groups = self._groups_widget()

        def update(_=None):
            self.points_filepath = wgt_file.value
            # the groups are grayed out as soon as a points file is typed
            self._set_enabled(wgt_groups, not self.points_filepath)
            wgt_file_info.value = self._points_file_info(self.points_filepath)

        wgt_file.observe(update, names="value")
        update()
        return widgets.VBox([
            self._section("Points file", wgt_file, wgt_file_info,
                          self._note("Columns: Group (optional), Distance (AU), Particle, Energy (MeV, e.g. 13.5 or "
                                     "13-16) and Onset Time, one row per point (see "
                                     "examples/custom_points_example.csv)")),
            self._section("Groups of points (used when there is no points file)", wgt_groups),
        ])

    def _points_file_info(self, filepath):
        if not filepath:
            return ""
        if not os.path.isfile(filepath):
            return "<i>File not found</i>"
        try:
            groups = read_points_file(filepath)
        except ValueError as e:
            return f"<pre>{html.escape(str(e))}</pre>"
        n_points = sum(len(vda.points) for vda in groups.values())
        return f"<i>{len(groups)} groups, {n_points} points: {html.escape(', '.join(groups))}</i>"

    def _groups_widget(self):
        """A tab per group of points, and a button to add groups"""
        tabs = widgets.Tab()

        def show_groups(selected=0):
            tabs.children = [self._group_widget(i, show_groups) for i in range(len(self.groups))]
            for i, group in enumerate(self.groups):
                tabs.set_title(i, group["label"])
            tabs.selected_index = selected

        def add_group(_):
            self.groups.append(self._new_group(len(self.groups) + 1))
            show_groups(len(self.groups) - 1)

        btn_add = widgets.Button(description="Add Group")
        btn_add.on_click(add_group)
        show_groups()
        self._groups_tabs = tabs
        return widgets.VBox([tabs, btn_add])

    def _group_widget(self, i, show_groups):
        """Label, distance and points of a group, with buttons to add and remove points"""
        group = self.groups[i]
        wrapper_rows = widgets.VBox()

        wgt_label = widgets.Text(value=group["label"], description="Label:", style=self._style)

        def set_label(traitlet):
            group["label"] = traitlet["new"]
            self._groups_tabs.set_title(i, traitlet["new"])

        wgt_label.observe(set_label, names="value")
        wgt_distance = widgets.BoundedFloatText(value=group["distance"], min=0.01, max=100, step=0.0001,
                                                description="Distance from the Sun (AU):", style=self._style)
        wgt_distance.observe(lambda traitlet: group.__setitem__("distance", traitlet["new"]), names="value")
        btn_remove_group = widgets.Button(description="Remove Group", tooltip=f"Remove {group['label']}",
                                          disabled=len(self.groups) == 1)

        def remove_group(_):
            del self.groups[i]
            show_groups(max(i - 1, 0))

        btn_remove_group.on_click(remove_group)

        def set_point_value(j, k, value):
            group["points"][j][k] = value

        def remove_point(j):
            del group["points"][j]
            show_rows()

        def add_point(_):
            # a copy of the last point, to be changed
            group["points"].append(list(group["points"][-1]) if group["points"] else ["", "", "protons"])
            show_rows()

        def show_rows():
            rows = []
            for j, (onset, energy, particle) in enumerate(group["points"]):
                wgt_onset = widgets.Text(value=str(onset), description=f"Point {j + 1} onset:", style=self._style,
                                         placeholder="e.g. 2024-11-21 18:20:14")
                wgt_onset.observe(lambda traitlet, j=j: set_point_value(j, 0, traitlet["new"]), names="value")
                wgt_energy = widgets.Text(value=str(energy), description="energy (MeV):", style=self._style,
                                          placeholder="e.g. 13.5 or 13-16")
                wgt_energy.observe(lambda traitlet, j=j: set_point_value(j, 1, traitlet["new"]), names="value")
                wgt_particle = widgets.Dropdown(options=["protons", "electrons"],
                                                value=PARTICLE_NAMES.get(str(particle).strip().lower(), "protons"),
                                                description="particle:", style=self._style)
                wgt_particle.observe(lambda traitlet, j=j: set_point_value(j, 2, traitlet["new"]), names="value")
                btn_remove = widgets.Button(description="Remove", tooltip=f"Remove point {j + 1}")
                btn_remove.on_click(lambda _, j=j: remove_point(j))
                rows.append(widgets.HBox([wgt_onset, wgt_energy, wgt_particle, btn_remove]))
            wrapper_rows.children = rows

        btn_add = widgets.Button(description="Add Point")
        btn_add.on_click(add_point)
        show_rows()
        return widgets.VBox([
            widgets.HBox([wgt_label, wgt_distance, btn_remove_group]),
            wrapper_rows,
            btn_add,
        ])

    def _plots_sections(self):
        wgt_title = widgets.Text(value=self.title, description="Title:", style=self._style, layout=self._layout,
                                 placeholder="Title of the VDA plot - Leave blank for 'VDA of' the date of the onsets")
        wgt_title.observe(lambda traitlet: setattr(self, "title", traitlet["new"]), names="value")
        plots, *sections = super()._plots_sections()
        plots.children += (wgt_title,)
        return [plots, *sections]

    def _default_reference_time(self):
        # the earliest onset of the groups
        onsets = pd.to_datetime([p[0] for g in self.groups for p in g["points"]], errors="coerce").dropna()
        return onsets.min().to_pydatetime().replace(second=0) if len(onsets) else super()._default_reference_time()

    ############### VDA ###############
    def group_vdas(self) -> dict:
        """The VDA fit of each group of points ({label: CustomPointsVDA}), from the points file or the form"""
        if self.points_filepath:
            return read_points_file(self.points_filepath)
        vdas = {}
        for group in self.groups:
            label = group["label"].strip()
            if not label or label in vdas:
                raise ValueError(f"Group {label!r}: give each group a different label")
            try:
                points = points_df(group["points"])
            except ValueError as e:
                raise ValueError(f"Group {label!r}: {e}") from None
            vdas[label] = CustomPointsVDA(points, group["distance"], name=f"group {label}")
        return vdas

    def run_vda(self):
        """Fits the VDA line of each group of points, prints the results and plots them together"""
        self.vdas = self.group_vdas()
        for vda in self.vdas.values():
            vda.compute_vda()
        self.results = pd.DataFrame([vda.results.loc[CustomPointsVDA.EVENT_NO] for vda in self.vdas.values()],
                                    index=pd.Index(list(self.vdas), name="Group"))
        for label, vda in self.vdas.items():
            _print_results(label, self.results.loc[label], vda.format_timedelta)
        filename = views.custom_points_plot_filename(self.vdas) if self.save_plots else None
        references = [(label, time) for label, time in self.reference_times if time is not None]
        _show_figure(views.plot_vda_combined(self.vdas, CustomPointsVDA.EVENT_NO, filename, references,
                                             self.title or None))
