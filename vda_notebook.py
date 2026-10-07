"""Notebook interface of the VDA tool: the parameters form and the steps of the analysis.

The analysis itself is in vda.py, its parameters in vda_tool_configuration.py, the plots in vda_views.py
and the widgets in vda_widgets.py.
"""
import html
import os
from warnings import filterwarnings, simplefilter

from IPython.display import display
from ipywidgets import widgets
from matplotlib import pyplot as plt
from pandas.errors import PerformanceWarning
from sunpy import log as sunpy_log

import vda_views
from vda import VDA
from vda_tool_configuration import OnsetSelection, VDA_parameters
from vda_widgets import VDA_nb_displayer, show_figure


class VDA_notebook:

    def __init__(self):
        simplefilter(action="ignore", category=PerformanceWarning)
        filterwarnings(action="ignore", message="Discarding nonzero nanoseconds in conversion")

        self.vda = VDA(VDA_parameters())
        self.displayer = VDA_nb_displayer(self.vda)
        self.font_size = 22
        self.legend_font_size = 18
        self.save_plots = True
        self._apply_fonts()

    @property
    def parameters(self):
        return self.vda.parameters

    def _apply_fonts(self):
        plt.rc("font", size=self.font_size)
        plt.rc("legend", fontsize=self.legend_font_size)

    def _show(self, df):
        if self.parameters.view_dfs:
            display(df)

    ############### Parameters ###############
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

    def _set_enabled(self, widget, enabled):
        """Grays out (or activates again) the widget and all the widgets in it"""
        if hasattr(widget, "disabled"):
            # a single event cannot be removed
            single_event = getattr(widget, "description", "") == "Remove Event" and len(self.parameters.date_ranges) == 1
            widget.disabled = not enabled or single_event
        for child in getattr(widget, "children", ()):
            self._set_enabled(child, enabled)

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
            return [self.displayer.display_reference_hours()]
        if kind == "error":
            return [widgets.HTML(f"<pre>{html.escape(text)}</pre>")]
        return [self._note(text)]

    def _events_tab(self):
        wgt_file = self.displayer.display_input_file()
        # the date ranges are grayed out as soon as an events file is typed
        wgt_file.continuous_update = True
        wgt_ranges = self.displayer.display_date_ranges()
        wrapper_file = widgets.VBox()
        shown = {"info": None}

        def update(_=None):
            filepath = self.parameters.input_filepath
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

    def _data_tab(self):
        return widgets.VBox([
            self._section("Viewings", self.displayer.display_viewings()),
            self._section("Resampling", self.displayer.display_resample_frequency()),
            self._section("Advanced: load previously saved data, or save the data",
                          self.displayer.display_load_data_option(),
                          self.displayer.display_save_data_option()),
        ])

    def _channels_tab(self):
        # the energy ranges of the channels are the same for all the events
        start, end = self.parameters.date_ranges[0]
        self.vda.construct_energies_df(start, end)
        return self._section(
            "Grouped energy channels (select multiple channels with Ctrl+click)",
            self.displayer.display_channel_groups(),
        )

    def _onsets_tab(self):
        return widgets.VBox([
            self._section("Onset determination method",
                          self.displayer.display_onset_method_selection(),
                          self.displayer.display_onset_method_parameters()),
            self._section("Default background window",
                          self.displayer.display_bg_defaults(),
                          self._note("The background window of each event can be checked and adjusted later, "
                                     "after the data are downloaded (Background windows step).")),
            self._section("Onset selection", self.displayer.display_onset_selection_selection()),
        ])

    def _plots_tab(self):
        wgt_font = widgets.BoundedIntText(value=self.font_size, min=6, max=40, description="Font size:",
                                          style=self.displayer.WIDGETS_STYLE)
        wgt_legend_font = widgets.BoundedIntText(value=self.legend_font_size, min=6, max=40, description="Legend font size:",
                                                 style=self.displayer.WIDGETS_STYLE)
        wgt_save = widgets.Checkbox(value=self.save_plots, description="Save the VDA plots as .png files", indent=False)

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
            self._section("Tables", self.displayer.display_view_toggle()),
        ])

    ############### Analysis ###############
    def run_data(self):
        """Creates the events, downloads (or loads) their data and groups the energy channels"""
        self.vda.construct_times_df()
        self._show(self.vda.df_times)
        self.vda.construct_particles_df()
        self.vda.group_energy_channels()
        self._show(self.vda.df_grouped)

    def background_selection(self):
        """Background window of each event, with a slider to change it"""
        return self.displayer.display_bg_selection()

    def run_onsets(self):
        """Determines the onsets and selects the ones used for the VDA (interactively, if selected)"""
        self.vda.calculate_onsets()
        self.vda.clean_onsets()
        self.vda.construct_options_df()
        self._show(self.vda.df_onsets_existing)
        self.vda.select_onsets()
        if self.parameters.onset_selection == OnsetSelection.INTERACTIVE:
            return self.displayer.display_onset_selection()

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
            filename = vda_views.vda_plot_filename(self.vda, event_no) if self.save_plots else None
            show_figure(vda_views.plot_vda(self.vda, event_no, filename))
        self._show(self.vda.results)
