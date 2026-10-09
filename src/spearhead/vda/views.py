"""Plots of the VDA analysis: background window, detected onsets and VDA fit of an event.

The functions return matplotlib (pyplot) figures without showing them, so they can be used in scripts, notebooks
or widgets. The figures are saved only when a filename is given, and can be closed with plt.close(fig) when they
are no longer needed.
"""
from datetime import timedelta
from math import ceil

import matplotlib
import numpy as np
from matplotlib import dates as mdates
from matplotlib import pyplot as plt
from matplotlib.figure import Figure


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
    bg_end_format = "%H:%M" if bg_start.date() == bg_end.date() else "%Y-%m-%d %H:%M"
    ax.set_title(
        f"Event {event_no}\n"
        f"Background: {bg_start:%Y-%m-%d %H:%M} to {bg_end:{bg_end_format}} ({vda.bg_window_source(event_no)})"
    )

    ax.set_yscale("log")
    ax.set_ylabel(r"Intensity (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$)")
    _time_axis(ax, temp_df.index)
    ax.set_xlabel("Time (UTC)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize="small")
    if filename:
        fig.savefig(filename)
    return fig


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
        values = [f"Onset {onset_results['Onset Time']:%H:%M}"]
        # values of the onset method: threshold (sigma), mu_d and the k and h parameters (Poisson-CUSUM)
        method_specific = onset_results["Method Specific"]
        if not isinstance(method_specific, dict):
            method_specific = {}
        if "bg_level" in method_specific:
            ax.axhline(method_specific["bg_level"], color="green", linestyle="dashed", label="BG level")
            values.append(f"BG {method_specific['bg_level']:.3g}")
        if "threshold" in method_specific:
            ax.axhline(method_specific["threshold"], color="red", linestyle="dashed", label="Threshold")
            values.append(f"Threshold {method_specific['threshold']:.3g}")
        if "mu_d" in method_specific:
            ax.axhline(method_specific["mu_d"], color="red", linestyle="dashed", label="μd")
            values += [f"μd {method_specific['mu_d']:.3g}",
                       f"k {method_specific['k']:.3g}, h {method_specific['h']:.3g}"]
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
        f"{sensor}/{particle} ({low_energy:.2f}-{high_energy:.2f} MeV)"
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


def plot_vda(vda, event_no, filename=None) -> Figure:
    """Plots the VDA fit of the event, computed by vda.compute_vda()"""
    fit = vda.vda_fits[event_no]
    inv_betas = fit["inv_betas"]
    a, b, b_error = fit["a"], fit["b"], fit["b_error"]

    def to_time(seconds):
        return fit["t0"] + timedelta(seconds=seconds)

    res = vda.results.loc[event_no]
    fig = plt.figure(figsize=(12, 8), layout="constrained")
    ax = fig.add_subplot()
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
    fig.suptitle(f"Event {event_no} ({vda.df_grouped.loc[event_no].index[1].to_pydatetime().strftime('%Y-%m-%d')})")
    ax.set_xlabel("Inverse Beta")
    ax.set_ylabel("Time")
    ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # the results are shown as legend entries without markers
    for label in (
        f"Extra Time = {vda.format_timedelta(res['Extra Time'])}",
        f"Release Time = {res['Release Time']} +/- {vda.format_timedelta(res['Release Time Error'])}",
        f"APL = {res['APL']:.2f} +/- {res['APL Error']:.2f}",
    ):
        ax.plot([], [], alpha=0, label=label)
    # legend between the title and the plot
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncols=2, frameon=False)
    if filename:
        fig.savefig(filename)
    return fig
