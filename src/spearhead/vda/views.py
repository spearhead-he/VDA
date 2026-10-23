"""Plots of the VDA analysis: background window, detected onsets and VDA fit of an event.

The functions return matplotlib (pyplot) figures without showing them, so they can be used in scripts, notebooks
or widgets. The figures are saved only when a filename is given, and can be closed with plt.close(fig) when they
are no longer needed.
"""
from datetime import timedelta
from itertools import cycle, islice
from math import ceil, sqrt

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import dates as mdates
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .analysis import format_error, format_timedelta, reference_time_at_sun, release_times, time_format


def _time_axis(ax, times, maxticks=None) -> None:
    """Concise datetime ticks for the data range of the event"""
    ax.set_xlim(times[0], times[-1])
    locator = mdates.AutoDateLocator() if maxticks is None else mdates.AutoDateLocator(maxticks=maxticks)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))


def plot_bg(vda, event_no, filename=None) -> Figure:
    """Plots the grouped channels of the event with its background window"""
    fig = plt.figure(figsize=(16, 7), layout="constrained")
    ax = fig.add_subplot()
    # one colormap per sensor and particle, darker for the later channels of the group
    colormaps = ["Blues", "Oranges", "Greens", "Purples", "Reds", "Greys"]
    linestyles = ["-", "--", ":", "-."]
    # zeros cannot be shown in log scale
    temp_df = vda.df_grouped.loc[event_no].replace(0, np.nan)
    for i_group, (sensor, particle, particle_prefix) in enumerate(vda._iter_sensor_particles()):
        cmap = matplotlib.colormaps[colormaps[i_group % len(colormaps)]]
        viewings = vda.sensor_viewings(sensor)
        for i_viewing, viewing in enumerate(viewings):
            df_channels = temp_df[sensor][particle][viewing][particle_prefix]
            n_channels = len(df_channels.columns)
            for i_channel, channel in enumerate(df_channels.columns):
                low_energy, high_energy = vda.channel_energy_range(sensor, particle, particle_prefix, channel)
                label = f"{sensor.upper()} {particle} {low_energy:.2f}-{high_energy:.2f} MeV"
                if len(viewings) > 1:
                    label += f" ({viewing})"
                ax.plot(
                    df_channels[channel],
                    color=cmap(0.4 + 0.6 * i_channel / max(n_channels - 1, 1)),
                    linestyle=linestyles[i_viewing % len(linestyles)],
                    label=label,
                )

    bg_start, bg_end = vda.bg_window(event_no)
    ax.axvspan(bg_start, bg_end, color="green", alpha=0.2, label="Background")
    bg_time_format = time_format(bg_start, bg_end)
    bg_end_format = bg_time_format if bg_start.date() == bg_end.date() else f"%Y-%m-%d {bg_time_format}"
    ax.set_title(
        f"Event {event_no}\n"
        f"Background: {bg_start:%Y-%m-%d {bg_time_format}} to {bg_end:{bg_end_format}} ({vda.bg_window_source(event_no)})"
    )

    ax.set_yscale("log")
    ax.set_ylabel(r"Intensity (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$)")
    _time_axis(ax, temp_df.index)
    ax.set_xlabel("Time (UTC)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize="small")
    if filename:
        fig.savefig(filename)
    return fig


# level of the onset methods, with its color and legend label: background level, threshold (sigma) and mu_d
# (Poisson-CUSUM)
ONSET_LEVELS = [("bg_level", "green", "BG level"), ("threshold", "red", "Threshold"), ("mu_d", "red", "μd")]


def _plot_onset_levels(ax, method_specific: dict) -> list[str]:
    """Plots the levels of the onset method, and returns the texts of its values (also the k and h parameters of
    Poisson-CUSUM)"""
    for key, color, label in ONSET_LEVELS:
        # a level of zero cannot be shown in log scale
        if method_specific.get(key, 0) > 0:
            ax.axhline(method_specific[key], color=color, linestyle="dashed", label=label)
    names = {"bg_level": "BG", "threshold": "Threshold", "mu_d": "μd"}
    values = [f"{names[key]} {method_specific[key]:.3g}" for key, _, _ in ONSET_LEVELS if key in method_specific]
    if "k" in method_specific:
        values.append(f"k {method_specific['k']:.3g}, h {method_specific['h']:.3g}")
    if method_specific.get("bg_level") == 0:
        values.append("Zero background (no counts):\nany counts are above it")
    return values


def plot_onsets(vda, event_no, sensor, particle, particle_prefix, channel, selected_viewing=None, filename=None) -> Figure:
    """Plots the detected onsets of a grouped channel of the event, one subplot per viewing of its sensor"""
    viewings = vda.sensor_viewings(sensor)
    temp_df = vda.df_grouped.loc[event_no]
    ncols = min(len(viewings), 3)
    nrows = ceil(len(viewings) / ncols)
    fig = plt.figure(figsize=(14, 4.5 * nrows + 1), dpi=300, layout="constrained")
    # the last row holds the legend of all the viewings
    grid = fig.add_gridspec(nrows + 1, ncols, height_ratios=[1] * nrows + [0.15])
    axs_flat = [fig.add_subplot(grid[row, col]) for row in range(nrows) for col in range(ncols)]
    ax_legend = fig.add_subplot(grid[nrows, :])
    ax_legend.axis("off")
    for ax in axs_flat[len(viewings):]:
        ax.axis("off")
    for ax, viewing in zip(axs_flat, viewings):
        ax.set_title(f"{viewing} (selected)" if viewing == selected_viewing else viewing,
                     fontweight="bold" if viewing == selected_viewing else "normal")
        try:
            onset_results = vda.df_onsets_existing.loc[(event_no, sensor, particle, viewing, particle_prefix, channel)]
        except KeyError:
            ax.text(0.5, 0.5, "No onset", transform=ax.transAxes, ha="center", va="center")
            continue

        ax.plot(temp_df[sensor][particle][viewing][particle_prefix][channel].fillna(0).ffill(), label="Data")
        ax.set_yscale("log")
        ax.axvspan(onset_results["Background Start"], onset_results["Background End"],
                   color="green", alpha=0.3, label="BG sample")
        method_specific = onset_results["Method Specific"]
        onset_time = onset_results["Onset Time"]
        values = [f"Onset {onset_time:{time_format(onset_time)}}"]
        values += _plot_onset_levels(ax, method_specific if isinstance(method_specific, dict) else {})
        ax.axvline(onset_results["Onset Time"], color="purple", linestyle="dashed", label="Onset")
        ax.text(0.98, 0.03, "\n".join(values), transform=ax.transAxes, ha="right", va="bottom",
                fontsize="x-small", bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
        _time_axis(ax, temp_df.index, maxticks=5)

    # one legend for all the viewings, from the first subplot with an onset
    for ax in axs_flat:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax_legend.legend(handles, labels, loc="center", ncols=len(handles), fontsize="small", frameon=False)
            break

    low_energy, high_energy = vda.channel_energy_range(sensor, particle, particle_prefix, channel)
    fig.suptitle(
        f"Detected onsets for event {event_no} ({temp_df.index[0]:%Y-%m-%d}) | "
        f"{sensor.upper()} {particle} ({low_energy:.2f}-{high_energy:.2f} MeV)"
    )
    if filename:
        fig.savefig(filename)
    return fig


def vda_plot_filename(vda, event_no) -> str:
    """Default filename of the VDA plot of the event: data range, sensors and particles, resample frequency"""
    event_times = vda.df_times.loc[event_no]
    time_start_str = event_times[vda.START_TIME_COLNAME].strftime("%Y-%m-%d_%H%M")
    time_end_str = event_times[vda.END_TIME_COLNAME].strftime("%Y-%m-%d_%H%M")
    particles_str = "_".join([f"{s}-{p}" for s, ps in vda.parameters.sensors_particles.items() for p in ps])
    freq_str = vda.parameters.resample_frequency if vda.parameters.resample_frequency != "" else "noresample"
    return f"{time_start_str}_{time_end_str}_{particles_str}_{freq_str}.png"


def _particle_markers() -> list:
    """(particle, marker, size) of the VDA points: protons as circles, electrons as stars, enlarged to look as big as
    the circles"""
    marker_size = plt.rcParams["lines.markersize"] ** 2
    return [("protons", "o", marker_size), ("electrons", "*", 2.25 * marker_size)]


# marker and color of each reference time, in order. The colors are not those of the spacecraft in the combined plot
REFERENCE_STYLES = [("^", "tab:pink"), ("v", "tab:olive"), ("X", "black"), ("D", "tab:cyan"), ("s", "tab:gray"),
                    ("P", "navy")]


def _plot_references(ax, references) -> None:
    """Reference times (e.g. SXR peak, type II and III radio bursts) at inverse beta 0, at the Sun, with one legend entry
    each giving the time at the Sun and the observed time, e.g. "SXR peak: 15:26:44 (15:35:00 at Earth)".

    references: (label, datetime) observed at Earth, or (label, datetime, observer), see
    analysis.reference_time_at_sun
    """
    marker_size = plt.rcParams["lines.markersize"] ** 2
    for reference, (marker, color) in zip(references, cycle(REFERENCE_STYLES)):
        label, observed, observer, at_sun = reference_time_at_sun(reference)
        ax.scatter([0], [at_sun.to_pydatetime(warn=False)], marker=marker, color=color, s=2 * marker_size,
                   label=f"{label}: {at_sun:%H:%M:%S} ({observed:%H:%M:%S} at {observer})", zorder=5, clip_on=False)
    ax.set_xlim(left=0)


def _plot_fit_band(ax,
                   fit_x: np.ndarray,
                   fit: dict,
                   color: str) -> bool:
    """Plots the ±2σ band of the release time around the fit line (in vda.vda_fits). Returns False if the fit has no
    error (2 points)"""
    a, b, b_error = fit["a"], fit["b"], fit["b_error"]
    if np.isnan(b_error):
        return False
    ax.fill_between(
        fit_x,
        [fit["t0"] + timedelta(seconds=a * x + b - 2 * b_error) for x in fit_x],
        [fit["t0"] + timedelta(seconds=a * x + b + 2 * b_error) for x in fit_x],
        color=color,
        alpha=0.1,
    )
    return True


def _vda_figure(n_fits: int) -> tuple:
    """Figure of a VDA plot: the plot, and below it the axes of the results table, with a row per fit"""
    fig = plt.figure(figsize=(12, 9 + 0.5 * n_fits), layout="constrained")
    grid = fig.add_gridspec(2, 1, height_ratios=[8, 0.6 + 0.55 * n_fits])
    ax = fig.add_subplot(grid[0])
    ax_table = fig.add_subplot(grid[1])
    ax_table.axis("off")
    ax.set_xlabel("Inverse Beta")
    ax.set_ylabel("Time (UTC)")
    ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    return fig, ax, ax_table


def _plot_fit(ax,
              fit: dict,
              max_inv_beta: float,
              color: str,
              points_color: str) -> tuple[set, bool]:
    """Plots the points of a VDA fit (in vda.vda_fits), protons as circles and electrons as stars, and its line and ±2σ
    band from inverse beta 0, where an open circle marks the release time at the Sun.

    Returns the particles of the points, and whether the band is plotted.
    """
    fit_x = np.array([0, max_inv_beta])

    def to_time(seconds):
        return fit["t0"] + timedelta(seconds=seconds)

    particles = set()
    for particle, marker, size in _particle_markers():
        is_particle = fit["particles"] == particle
        if is_particle.any():
            particles.add(particle)
            ax.scatter(fit["inv_betas"][is_particle],
                       [to_time(t) for t in fit["onset_seconds"][is_particle]],
                       color=points_color, marker=marker, s=size, zorder=2)
    ax.plot(fit_x, [to_time(fit["a"] * x + fit["b"]) for x in fit_x], color=color)
    has_band = _plot_fit_band(ax, fit_x, fit, color)
    ax.scatter([0], [to_time(fit["b"])], marker="o", facecolor="white", edgecolor=color,
               linewidths=2, s=1.3 * plt.rcParams["lines.markersize"] ** 2, zorder=4, clip_on=False)
    return particles, has_band


# size of the legends inside the VDA plots, relative to the legend font size, to hide less of the plot
INNER_LEGEND_SCALE = 0.9


def _inner_legend_size() -> float:
    """Font size of the legends inside the VDA plots, in points"""
    return INNER_LEGEND_SCALE * FontProperties(size=plt.rcParams["legend.fontsize"]).get_size_in_points()


def _plot_legends(ax,
                  particles: set,
                  has_band: bool,
                  n_fits: int,
                  color: str) -> None:
    """Legends inside the plot: the symbols in the upper left corner (the points, the fit lines, their bands and the
    release times at the Sun), and the reference times (the labelled points of the plot) in the lower right corner"""
    plural = "s" if n_fits > 1 else ""
    handles = [Line2D([], [], color="black", marker=marker, linestyle="", markersize=sqrt(size),
                      label=particle.capitalize())
               for particle, marker, size in _particle_markers() if particle in particles]
    handles.append(Line2D([], [], color=color, label=f"VDA fit{plural}"))
    if has_band:
        handles.append(Patch(color=color, alpha=0.2, label=f"2σ band{plural}"))
    handles.append(Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=color, markeredgewidth=2,
                          linestyle="", markersize=1.15 * plt.rcParams["lines.markersize"], label="Release at the Sun"))
    references = ax.get_legend_handles_labels()[0]
    if references:
        ax.add_artist(ax.legend(handles=references, loc="lower right", fontsize=_inner_legend_size()))
    ax.legend(handles=handles, loc="upper left", fontsize=_inner_legend_size())


# columns of the results table under the VDA plots
RESULTS_COLUMNS = ["", "Release at the Sun", "At 1 AU", "Error", "APL (AU)"]


def results_rows(fits: list[tuple]) -> list[list[str]]:
    """Rows of the results table of a VDA plot, from (label, fit, results) of each fit (a fit of vda.vda_fits and its
    row of vda.results): the release time at the Sun, the time its light reaches 1 AU, their error (1σ) and the APL"""
    rows = []
    for label, fit, res in fits:
        at_sun, at_au = release_times(fit)
        rows.append([label, f"{at_sun:%H:%M:%S}", f"{at_au:%H:%M:%S}",
                     f"± {format_timedelta(res['Release Time Error'])}",
                     f"{res['APL']:.2f} ± {format_error(res['APL Error'])}"])
    return rows


def _results_table(ax_table,
                   rows: list[list[str]],
                   colors: list[str]) -> None:
    """Table of the results under the plot, with the label of each fit in its color"""
    if not rows:
        return
    table = ax_table.table(cellText=rows, colLabels=RESULTS_COLUMNS, loc="upper center", cellLoc="center",
                           colWidths=[0.22, 0.23, 0.15, 0.17, 0.21], edges="horizontal")
    table.auto_set_font_size(False)
    table.set_fontsize(plt.rcParams["legend.fontsize"])
    table.scale(1, 2)
    for row, color in enumerate(colors, start=1):
        table[row, 0].get_text().set_color(color)
        table[row, 0].get_text().set_fontweight("bold")


def plot_vda(vda, event_no, filename=None, references=None) -> Figure:
    """Plots the VDA fit of the event, computed by vda.compute_vda().

    The fit line and its ±2σ band are extended to inverse beta 0, where an open circle marks the release time at the
    Sun. The table under the plot gives it, the time its light reaches 1 AU, their error and the APL.
    references: optional list of reference times (e.g. SXR peak, type II and III radio bursts), (label, datetime)
    observed at Earth or (label, datetime, observer), with the observer "earth", a name of OBSERVERS (e.g. "sta") or
    its distance from the Sun in AU. They are shown at inverse beta 0 at the Sun, i.e. minus the light travel time
    from the Sun to their observer, with this time and the observed one in the legend
    """
    fit = vda.vda_fits[event_no]
    fig, ax, ax_table = _vda_figure(1)
    particles, has_band = _plot_fit(ax, fit, max(fit["inv_betas"]), "blue", "black")
    if references:
        _plot_references(ax, references)
    _plot_legends(ax, particles, has_band, 1, "blue")
    ax.set_xlim(0, 1.03 * max(fit["inv_betas"]))
    fig.suptitle(f"Event {event_no} ({vda.df_grouped.loc[event_no].index[1].to_pydatetime().strftime('%Y-%m-%d')})"
                 f" | {vda.observer.label}")
    _results_table(ax_table, results_rows([(vda.observer.label, fit, vda.results.loc[event_no])]), ["black"])
    if filename:
        fig.savefig(filename)
    return fig


def vda_combined_plot_filename(vdas, event_no=1) -> str:
    """Default filename of the combined VDA plot of the event: data range and observers"""
    vda = next(iter(vdas.values()))
    event_times = vda.df_times.loc[event_no]
    time_start_str = event_times[vda.START_TIME_COLNAME].strftime("%Y-%m-%d_%H%M")
    time_end_str = event_times[vda.END_TIME_COLNAME].strftime("%Y-%m-%d_%H%M")
    observers_str = "-".join(vda.parameters.observer for vda in vdas.values())
    return f"{time_start_str}_{time_end_str}_{observers_str}.png"


def plot_vda_combined(vdas, event_no=1, filename=None, references=None, title=None) -> Figure:
    """Plots together the VDA fits of the same event from several spacecraft, computed by vda.compute_vda() of each.

    vdas: {label: VDA}, e.g. {"STEREO-A": vda_sta, "SOHO": vda_soho}. The spacecraft without a fit of the event are
    left out. Each spacecraft has its color, with the proton points as circles and the electron points as stars, and
    its fit line and ±2σ band extended to inverse beta 0, where an open circle marks its release time at the Sun. The
    table under the plot gives, for each spacecraft, the release time at the Sun, the time its light reaches 1 AU,
    their error and the APL.
    references: optional list of reference times, as in plot_vda
    title: title of the plot, by default "VDA of" the date of the event
    """
    fits = {label: vda for label, vda in vdas.items() if event_no in getattr(vda, "vda_fits", {})}
    fig, ax, ax_table = _vda_figure(len(fits))
    max_inv_beta = max((max(vda.vda_fits[event_no]["inv_betas"]) for vda in fits.values()), default=1)
    colors = list(islice(cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"]), len(fits)))
    particles, has_band = set(), False
    for vda, color in zip(fits.values(), colors):
        fit_particles, fit_band = _plot_fit(ax, vda.vda_fits[event_no], max_inv_beta, color, color)
        particles |= fit_particles
        has_band |= fit_band
    if references:
        _plot_references(ax, references)
    if fits:
        _plot_legends(ax, particles, has_band, len(fits), "gray")
    ax.set_xlim(0, 1.03 * max_inv_beta)
    if title:
        fig.suptitle(title)
    elif fits:
        fig.suptitle(f"VDA of {next(iter(fits.values())).vda_fits[event_no]['t0']:%Y-%m-%d}")
    rows = results_rows([(label, vda.vda_fits[event_no], vda.results.loc[event_no]) for label, vda in fits.items()])
    _results_table(ax_table, rows, colors)
    if filename:
        fig.savefig(filename)
    return fig


def custom_points_plot_filename(vdas) -> str:
    """Default filename of the VDA plot of custom points (CustomPointsVDA of each group): the earliest onset"""
    onsets = [onset for vda in vdas.values() for onset in pd.to_datetime(vda.points["Onset Time"])]
    return f"{min(onsets):%Y-%m-%d_%H%M}_custom_points.png" if onsets else "custom_points.png"
