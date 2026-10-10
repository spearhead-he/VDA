# Using the VDA tool without the notebook

The analysis code can be used in scripts or other notebooks:

| File | Contents |
|---|---|
| `spearhead.vda.analysis` | `VDA`: the analysis (events, data, channel grouping, onsets, VDA fit) |
| `spearhead.vda.conf` | `VDA_parameters`: the parameters of the analysis, and the `OnsetSelection` enum |
| `spearhead.vda.observers` | The observers (spacecraft): their sensors, channels and viewings, data loading and distance from the Sun |
| `spearhead.vda.views` | The plots, as matplotlib figures (not displayed, saved only if a filename is given) |
| `spearhead.vda.notebook` | `VDA_notebook`: the interface of `vda_tool.ipynb` (parameters form, analysis steps and their widgets) |

`VDA`, `VDA_parameters` and `OnsetSelection` can also be imported directly from `spearhead.vda`.

## Installation

From the folder of the tool (see the [README](../README.md#install-locally)):

```bash
pip install -e .
```

or directly from GitHub, without the notebook:

```bash
pip install "spearhead-vda @ git+https://github.com/spearhead-he/VDA"
```

A specific release is installed with its tag, e.g. `pip install "spearhead-vda @ git+https://github.com/spearhead-he/VDA@v0.6.0"`.

The extra `[notebook]` (e.g. `pip install -e ".[notebook]"`) also installs the notebook dependencies (Jupyter and ipywidgets), needed by `spearhead.vda.notebook`. The extra `[test]` installs the tools to run the tests (`pytest tests/`). The installed version is `spearhead.vda.__version__`.

## Example

```python
from spearhead.vda import VDA, VDA_parameters, views

parameters = VDA_parameters()
parameters.input_filepath = "examples/datetime_range_example.csv"
parameters.observer = "solo"            # Solar Orbiter (or "sta", "psp", "soho", "wind", "bepi")
parameters.viewings = ["sun"]
parameters.channel_groups = {
    "protons": {
        "HET protons 1": {"sensor": "het", "channels": [1, 2, 3]},
        "HET protons 2": {"sensor": "het", "channels": [10, 11, 12]},
        "HET protons 3": {"sensor": "het", "channels": [19, 20, 21]},
    },
    "electrons": {
        "HET electrons 1": {"sensor": "het", "channels": [0, 1]},
        "HET electrons 2": {"sensor": "het", "channels": [2, 3]},
    },
}
parameters.view_dfs = False

vda = VDA(parameters)
vda.construct_times_df()                # events and background windows
vda.construct_energies_df()             # energy bins of the channels
vda.construct_particles_df()            # download (or load) the data
vda.group_energy_channels()
vda.set_bg_window(2, "2021-11-09 15:30", "2021-11-09 16:30")   # optional: background of one event
vda.calculate_onsets()
vda.clean_onsets()
vda.construct_options_df()
vda.select_onsets()                     # "Use all" onset selection
vda.construct_energy_channels_characteristics()
vda.define_spacecraft_parameters()      # distance of the observer from the Sun (e.g. SPICE kernels for Solar Orbiter)
vda.compute_vda()                       # VDA fit of every event, without plots
vda.print_results()
print(vda.results)

fig = views.plot_vda(vda, 1, filename="event_1.png")   # VDA plot of event 1
```

## Plots

The functions of `spearhead.vda.views` return a matplotlib figure without showing it, saved to `filename` if one is given:

| Function | Plot |
|---|---|
| `plot_bg(vda, event_no)` | Grouped channels of the event with its background window |
| `plot_onsets(vda, event_no, sensor, particle, prefix, channel, selected_viewing=None)` | Detected onsets of a grouped channel, one subplot per viewing |
| `plot_vda(vda, event_no)` | VDA fit of the event (after `vda.compute_vda()`) |

`vda_plot_filename(vda, event_no)` gives the default filename of the VDA plot. In a notebook, the returned figure is shown by Jupyter; in a script, it is saved with `fig.savefig(...)` or shown with `plt.show()`. Figures that are no longer needed are closed with `plt.close(fig)`.

## Widgets in another notebook

`VDA_notebook` can use an existing `VDA` object, to check the background windows or select the onsets interactively in another notebook:

```python
from spearhead.vda.notebook import VDA_notebook

tool = VDA_notebook(vda)
tool.background_selection()   # after vda.group_energy_channels()
tool.onset_selection()        # after vda.select_onsets()
```

## Parameters

`VDA_parameters` is a dataclass: `print(parameters)` shows all the parameters, and `dataclasses.replace(parameters, input_filepath="other.csv")` creates a modified copy. The main ones:

| Parameter | Meaning |
|---|---|
| `input_filepath` | Events file (see the [user guide](user_guide.md#events)). If empty, `date_ranges` is used |
| `date_ranges` | List of (start, end) datetimes, one per event, used without an events file |
| `bg_hours_prior`, `bg_hours_after` | Data range of the files with reference times, in hours before and after the reference time |
| `bg_after_start` | Default background window, (start, end) in minutes after the start of the data of each event |
| `observer` | Spacecraft of the data, one of `observers.OBSERVERS`: `"solo"` (Solar Orbiter EPD: HET and EPT), `"sta"` (STEREO-A IMPACT: HET and SEPT) `"psp"` (Parker Solar Probe ISOIS: EPI-Hi HET), `"soho"` (SOHO: ERNE-HED and COSTEP-EPHIN), `"wind"` (Wind: 3DP) or `"bepi"` (BepiColombo: SIXS-P), see the [user guide](user_guide.md#observers). Setting `observer` sets `viewings` and `channel_groups` to its defaults |
| `viewings` | Viewings, in the priority order of the "Use all" onset selection. Each sensor uses the selected viewings it has (`vda.sensor_viewings(sensor)`). By default, those of the observer (`["sun"]` for Solar Orbiter, `["sun", "omni"]` for STEREO-A, `["A"]` for Parker Solar Probe, `["omni"]` for SOHO and Wind, `["side0"]` for BepiColombo) |
| `resample_frequency` | Pandas offset alias (e.g. `"5min"`), or `""` for no resampling |
| `channel_groups` | Grouped channels per particle: `{label: {"sensor": e.g. "het", "channels": [...]}}`. By default, those of the observer (for Solar Orbiter, 8 groups of HET protons and 2 of HET electrons; for STEREO-A, 4 of HET protons and 1 of HET electrons; for Parker Solar Probe, 3 of HET protons and 9 of HET electrons; for SOHO, 2 of ERNE protons and 1 of EPHIN electrons; for Wind, 3 of 3DP protons and 3 of 3DP electrons; for BepiColombo, 3 of SIXS protons and 3 of SIXS electrons) |
| `onset_method`, `onset_method_parameters` | Onset determination method, `"sigma"` (parameters `s`, `n`) or `"poisson_cusum"` (parameters `cusum_minutes`, `sigma_multiplier`), see the [user guide](user_guide.md#onset-determination). Setting `onset_method` sets its parameters to their defaults |
| `onset_selection` | `OnsetSelection.USE_ALL` or `OnsetSelection.INTERACTIVE` |
| `load_data_filepath`, `save_data_filepath` | .pkl files to load the data from instead of downloading them, and to save them to. Empty means not used |

Setting an unknown parameter (e.g. a typo) or a removed one raises an error.

## Background windows

- `vda.set_bg_window(event_no, start, end)`: background window of one event
- `vda.reset_bg_window(event_no)`: back to the default window
- `vda.set_default_bg_window(start_minutes, end_minutes)`: default window, for the events without one from the events file or `set_bg_window`
- `vda.save_times("path.csv")`: saves the data range and background window of every event as an events file
