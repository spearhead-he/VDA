# VDA tool user guide

- [Notebook](#notebook)
- [Parameters](#parameters)
- [Events](#events)
- [Background window](#background-window)
- [Onset determination](#onset-determination)
- [Onset selection](#onset-selection)
- [Results](#results)

## Notebook

The notebook `vda_tool.ipynb` has a setup cell, a form with all the parameters, and one cell per step of the analysis:

```python
from spearhead.vda.notebook import VDA_notebook
tool = VDA_notebook()        # setup

tool.parameters_form()       # all the parameters
tool.run_data()              # events, download (or load) of the data, grouping of the energy channels
tool.background_selection()  # background window of each event
tool.run_onsets()            # onsets, and their selection
tool.run_vda()               # results and plots
```

Run the cells in order. If a parameter is changed, rerun the steps after it. For advanced use, `tool.vda` is the analysis object and `tool.parameters` its parameters (see [Using the VDA tool without the notebook](library.md)).

## Observers

The data of one spacecraft (observer) are analysed at a time:

| Observer | Sensors | Viewings | Distance from the Sun |
|---|---|---|---|
| Solar Orbiter | EPD/HET and EPD/EPT (protons, electrons) | sun, asun, north, south, omni | SPICE kernels |
| STEREO-A | IMPACT/HET (protons, electrons) and IMPACT/SEPT (ions, shown as protons, and electrons) | HET: omni; SEPT: sun, asun, north, south | JPL Horizons |

Each sensor uses the selected viewings it has: with STEREO-A, HET uses omni and SEPT the others, so the grouped channels of a sensor need at least one of its viewings selected. The distance from the Sun gives the light travel time of the VDA (Extra Time); JPL Horizons needs an internet connection.

The default grouped channels of STEREO-A are its consecutive HET channels, in groups of three for protons and of two for electrons; a single channel left at the end joins the previous group.

## Parameters

The parameters form has five tabs:

- **Events**: the events file, or the datetime ranges of the events if no file is given (see [Events](#events)). For files with reference times, the data range is set with the hours prior to and after the reference time.
- **Data**:
  - the observer (spacecraft). Changing it sets the viewings and the grouped channels to the defaults of the observer
  - the viewings
  - the resample frequency ([offset alias](https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases), blank for no resampling)
  - optionally, a .pkl file to load previously saved data from (instead of downloading them) and a .pkl file to save the data to. The resample frequency is saved with the data and used when they are loaded.
- **Energy channels**: the grouped channels ("Add Channel" / "Remove Channel"; select multiple channels with Ctrl+click). The lists show the energy range of each channel. The intensity of a group is the bin-width weighted mean $I'=\frac{\sum I_n \Delta E_n}{\sum \Delta E_n}$.
- **Onsets**: the onset determination method and its parameters (see [Onset determination](#onset-determination)), the default background window, and the onset selection method.
- **Views/Plots**: font sizes, saving of the VDA plots as .png files, and display of the intermediate tables.

## Events

The events are given with a .csv file, or as datetime ranges entered in the form ("Add Event" button) if no file is given. The first column of the file is the event number, and the type of the file is deduced from the names of the other columns:

| Columns after the event number | Meaning | Example |
|---|---|---|
| `Start Time, End Time` | Datetime range of the data. The default background window is used | [datetime_range_only_example.csv](../examples/datetime_range_only_example.csv) |
| `Start Time, BG End, End Time` | As above, with the background window from `Start Time` to `BG End` | [datetime_range_example.csv](../examples/datetime_range_example.csv) |
| `Start Time, BG Start, BG End, End Time` | As above, with the background window from `BG Start` to `BG End`. Empty cells use the default background window | [datetime_range_bg_example.csv](../examples/datetime_range_bg_example.csv) |
| `Reference Time` | Reference datetime of the event. The datetime range of the data is set with the hours prior to and after it | [reference_times_example.csv](../examples/reference_times_example.csv) |

The datetimes can be in any format supported by `pandas.to_datetime` (e.g. `2024-12-31 00:00:00`).

## Background window

Each event has a background window, used for the onset determination:
- the background window of the events file, if given,
- otherwise the default background window, in minutes after the start of the data of the event (Onsets tab).

After the data are downloaded, the background window of each event can be checked and changed with the event dropdown and the "Background" slider, or set back to the default one ("Reset to default"). The chosen windows can be saved with `tool.vda.save_times("path.csv")`, and the saved file can be used as the events file of later runs.

A warning is shown for background windows outside the data of the event, or with fewer than 3 data points.

## Onset determination

The onset of each grouped channel and viewing is determined with the background window of its event, with one of the methods:

- **Sigma threshold** (`sigma`): the threshold is the background mean plus `s` standard deviations of the background, and the onset is found when `n` consecutive points are above it.
- **Poisson-CUSUM** (`poisson_cusum`): the Poisson-CUSUM method of [pyonset](https://github.com/Christian-Palmroos/PyOnset). The intensities are standardized with the background mean and standard deviation, and the onset is searched after the background window: it is found when the CUSUM function stays above its threshold for `cusum_minutes` (converted to data points with the cadence of the data). `sigma_multiplier` sets μd, the background mean plus `sigma_multiplier` standard deviations, used for the k parameter of the CUSUM function.

With both methods, the onset is the last point before the rise: before the first of the `n` points above the threshold (Sigma threshold), or before the CUSUM function first exceeds its threshold (Poisson-CUSUM, as in pyonset).

The plots of the interactive onset selection show the background level and the threshold (Sigma threshold) or μd (Poisson-CUSUM), and their values. Changing the method sets its parameters to their defaults.

## Onset selection

- **Use all**: for each grouped energy channel, the onset of the first viewing (in the order of the viewings) with a determined onset is used. With the checkboxes the order is sun, asun, north, south, omni; a different priority can be set with `tool.parameters.viewings`, e.g. `["north", "sun"]`.
- **Interactive**: the determined onsets are shown one channel at a time, chosen with the event and channel dropdowns or the "Previous" / "Next" buttons, and the viewing whose onset is used is selected per channel (or "None", to leave the channel out). The selection starts from the viewings of "Use all".

## Results

For each event, the release time, the extra time (light travel time from the Sun to the spacecraft) and the apparent path length (APL) are printed, followed by the VDA plot, which can also be saved as a .png file (Views/Plots tab). The results of all the events are stored in the `tool.vda.results` table. For many events, `tool.vda.compute_vda()` followed by `tool.vda.print_results()` gives the results without creating the plots.
