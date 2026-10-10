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
| Parker Solar Probe | ISOIS/EPI-Hi HET (protons, electrons) | A (sunward, as sun), B (anti-sunward, as asun) | JPL Horizons |
| SOHO | ERNE-HED (protons) and COSTEP-EPHIN (electrons) | omni | JPL Horizons |
| Wind | 3DP (protons, electrons) | omni | JPL Horizons |
| BepiColombo | SIXS-P (protons, electrons) | side0 (perpendicular to the Sun direction), side1 and side2 (anti-sunward) | JPL Horizons |

Each sensor uses the selected viewings it has: with STEREO-A, HET uses omni and SEPT the others, so the grouped channels of a sensor need at least one of its viewings selected. The distance from the Sun gives the light travel time of the VDA (Extra Time); JPL Horizons needs an internet connection.

The default grouped channels of STEREO-A, Parker Solar Probe, SOHO, Wind and BepiColombo are their consecutive channels, in groups of three for protons and of two for electrons; a single channel left at the end joins the previous group.

Parker Solar Probe:
- the data are the 1-minute EPI-Hi HET rates (`PSP_ISOIS-EPIHI_L2-HET-RATES60`), available from mid-2021
- the proton channels are 3 to 11 (11.3-53.8 MeV): the other channels of the data are empty
- the electron data are count rates (the electron fluxes of the data have gaps during events), divided by the width of the channel (counts s⁻¹ MeV⁻¹). They are proportional to the intensity, but their values are not intensities

SOHO:
- the protons are the 1-minute ERNE-HED intensities (`SOHO_ERNE-HED_L2-1MIN`), channels 0 to 6 (13-64 MeV). Its three highest channels are not corrected and are left out, as in seppy
- the electrons are the 1-minute EPHIN level 2 intensities (from the University of Kiel), channels 0 (E150, 0.25-0.7 MeV) and 2 (E1300). The energy range of E1300 depends on the date: 2.64-10.4 MeV before the failure mode D of EPHIN (4 October 2017), 0.67-10.4 MeV since. E300 (channel 1) is deactivated since then and E3000 (channel 3) has no data, so both are left out
- the energy ranges of the channels are read from the first event, so an analysis with EPHIN channels cannot have events on both sides of 4 October 2017 (it raises an error): the events before and after it are analysed separately

Wind:
- the data are the omnidirectional fluxes of the 3DP solid state telescopes: protons from SST Open (`WI_SOSP_3DP`, 9 channels, about 70 keV-6.8 MeV) and electrons from SST Foil (`WI_SFSP_3DP`, 7 channels, about 27-520 keV), at about 12 s
- the energy of each channel is its mean energy in the loaded data, with a width of 60% of it, as in seppy: the energy ranges of neighbouring channels overlap

BepiColombo:
- the data are the SIXS-P level 3 cruise phase data of the SERPENTINE project (monthly files of about 40 MB, at 2 minutes), available from 2020 to April 2024 with gaps: protons P1-P9 (about 1-90 MeV) and electrons E1-E7 (about 55 keV-10 MeV)
- the viewings are the detector sides. In the cruise phase the Sun is along the +Y axis of the spacecraft: about 97° from the boresight of side 0, and 135° from sides 1 and 2. Side 3, at 45° from the Sun, has no data in the level 3 product, and side 4 is blocked by the sunshade, so there is no sunward viewing
- the energy ranges of the channels differ slightly between the sides (most for P9 and E7); those of side 0 are used for all the sides

## Parameters

The parameters form has five tabs:

- **Events**: the events file, or the datetime ranges of the events if no file is given (see [Events](#events)). For files with reference times, the data range is set with the hours prior to and after the reference time.
- **Data**:
  - the observer (spacecraft). Changing it sets the viewings and the grouped channels to the defaults of the observer. The energy ranges of its channels are read from the data of the first event: if there are no data (e.g. BepiColombo before 2020), an error is shown in red and the grouped channels (Energy channels tab) are grayed out until another spacecraft is chosen, or the events are changed and the channels are read again ("Read the channels again" button); the rest of the form can still be used
  - the viewings
  - the resample frequency ([offset alias](https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases), blank for no resampling)
  - optionally, a .pkl file to load previously saved data from (instead of downloading them) and a .pkl file to save the data to. The resample frequency is saved with the data and used when they are loaded.
- **Energy channels**: the grouped channels ("Add Channel" / "Remove Channel"; select multiple channels with Ctrl+click). The lists show the energy range of each channel. The intensity of a group is the bin-width weighted mean $I'=\frac{\sum I_n \Delta E_n}{\sum \Delta E_n}$. A note shows when the energy ranges of a sensor changed on a date (SOHO EPHIN, see [Observers](#observers)).
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
