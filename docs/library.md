# Using the VDA tool without the notebook

The analysis code can be used in scripts or other notebooks:

| File | Contents |
|---|---|
| `vda.py` | `VDA`: the analysis (events, data, channel grouping, onsets, VDA fit, plots) |
| `vda_tool_configuration.py` | `VDA_parameters`: the parameters of the analysis, and the `OnsetSelection` enum |
| `vda_views.py` | `VDA_nb_displayer`: the widgets and the onset selection |
| `vda_notebook.py` | `VDA_notebook`: the interface of `vda_tool.ipynb` (parameters form and analysis steps) |

The files are not an installable package yet: run the scripts from the folder of the tool, or add it to the Python path.

## Example

```python
from vda import VDA
from vda_tool_configuration import VDA_parameters
from vda_views import VDA_nb_displayer

parameters = VDA_parameters()
parameters.input_filepath = "examples/datetime_range_example.csv"
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
VDA_nb_displayer(vda).select_onsets()   # "Use all" onset selection
vda.construct_energy_channels_characteristics()
vda.define_spacecraft_parameters()      # SPICE kernels for the spacecraft distance
vda.compute_vda()                       # VDA fit of every event, without plots
vda.print_results()
print(vda.results)
```

`vda.plot()` instead of `vda.compute_vda()` also creates the VDA plots.

## Parameters

`VDA_parameters` is a dataclass: `print(parameters)` shows all the parameters, and `dataclasses.replace(parameters, input_filepath="other.csv")` creates a modified copy. The main ones:

| Parameter | Meaning |
|---|---|
| `input_filepath` | Events file (see the [user guide](user_guide.md#events)). If empty, `date_ranges` is used |
| `date_ranges` | List of (start, end) datetimes, one per event, used without an events file |
| `bg_hours_prior`, `bg_hours_after` | Data range of the files with reference times, in hours before and after the reference time |
| `bg_after_start` | Default background window, (start, end) in minutes after the start of the data of each event |
| `viewings` | Viewings, in the priority order of the "Use all" onset selection |
| `resample_frequency` | Pandas offset alias (e.g. `"5min"`), or `""` for no resampling |
| `channel_groups` | Grouped channels per particle: `{label: {"sensor": "het" or "ept", "channels": [...]}}` |
| `onset_method`, `onset_method_parameters` | Onset determination method (`"sigma"`) and its parameters (`s`, `n`) |
| `onset_selection` | `OnsetSelection.USE_ALL` or `OnsetSelection.INTERACTIVE` |
| `load_data_filepath`, `save_data_filepath` | .pkl files to load the data from instead of downloading them, and to save them to. Empty means not used |

Setting an unknown parameter (e.g. a typo) or a removed one raises an error.

## Background windows

- `vda.set_bg_window(event_no, start, end)`: background window of one event
- `vda.reset_bg_window(event_no)`: back to the default window
- `vda.set_default_bg_window(start_minutes, end_minutes)`: default window, for the events without one from the events file or `set_bg_window`
- `vda.save_times("path.csv")`: saves the data range and background window of every event as an events file
