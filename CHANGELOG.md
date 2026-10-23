# Changelog

All notable changes of the SPEARHEAD VDA tool. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- The release time at the Sun and the time its light reaches 1 AU, as in the VDA plots, in the printed results and in the results table (columns `Release Time (Sun)` and `Release Time (1 AU)`). The printed results give the release time at the Sun with its error, then the times at 1 AU and at the spacecraft (or observer), with the light travel time from the Sun that the latter includes, e.g. `At the spacecraft : 2021-10-28 15:28:41 (+0:06:40 light travel time from the Sun)`. The user guide explains the release times, their error and the band of the plot (see [Results](docs/user_guide.md#results)).
- `vda.run()`: runs the steps of the analysis in scripts, from the events to the VDA fits, and returns the results. The tables of the steps stay available (e.g. `vda.df_onsets_existing`), and `bg_windows` sets the background window of some events, e.g. `vda.run(bg_windows={2: ("2021-10-28 14:10", "2021-10-28 14:55")})`. The steps can still be run one by one (see [Using the VDA tool without the notebook](docs/library.md#example)).
- `analysis.release_times(fit)`: the release time at the Sun of a VDA fit, where its line meets inverse beta 0, and the time its light reaches 1 AU.
- The observer of each reference time: Earth (the default), one of the spacecraft of the tool, or its distance from the Sun in AU, which is used instead of the observer ("Observed from" and "or distance (AU)" in the Views/Plots tab; `(label, datetime, observer)` in scripts, e.g. `("type III", datetime(2021, 10, 28, 15, 27), "sta")` or `0.8`). The reference times are shown at the Sun, minus the light travel time from the Sun to their observer at that time, and their legend, in the lower right corner of the plot, gives the time at the Sun and the observed time of each, e.g. `SXR peak: 15:18:44 (15:27:00 at Earth)` (`analysis.reference_time_at_sun`).
- Warnings with the results of a VDA fit that is not physical: an APL that is not positive, or a release at the Sun after the earliest onset (`analysis.fit_warnings`, in `vda.vda_fits[event_no]["warnings"]`).
- The Onsets step prints the number of grouped channels of each event with an onset, e.g. `Event 1: onsets in 8 of 10 grouped channels`, or that an event has no onsets (`vda.onset_counts()`, `vda.print_onset_counts()`).
- Clear error messages, before the data are downloaded where possible, for: a step run before the previous one (e.g. "Run the data step (tool.run_data()) first"), an events, points or data file that is not found, invalid datetimes in an events file (with their rows), an events file with a start after the end, a background window set for an event that does not exist, a grouped channel without channels or no grouped channels, an invalid resample frequency (also shown in red in the form), a points or events file whose columns are separated by `;` or tabs instead of commas, and a distance of an observer (group of points, reference time) that is not above 0 and at most 10 AU (e.g. a distance in km). The observer of a reference time can also be given by its label, in any case (e.g. `"STEREO-A"`).
- The parameters form says when it reads the energy ranges of the channels from the data of the first event, which may be downloaded.

### Changed
- Resampled data are timestamped at the centre of their bins, as the data loaders (seppy, solo-epd-loader) timestamp the data, instead of at their start (rounded down to the minute). The bins start at the start time of the event, for all the sensors and viewings (`analysis.resample`). The onset times, and so the release times, are half a bin later (2.5 minutes for the default 5 minutes), and a background window from 14:00 to 15:00 has the 12 bins inside it instead of 13. For the default event (2021-10-28), the release times at the spacecraft are Solar Orbiter 15:28:41 (instead of 15:26:11), STEREO-A 15:30:53 (15:28:23) and SOHO 15:44:00 ± 0:13:18 with an APL of 1.54 ± 0.44 AU (15:39:34 ± 0:05:14, 1.43 ± 0.17 AU): the onsets of two SOHO channels, close to their threshold, are one and two bins later with the bins and the background window. The background window slider and status, the title of the background plot and the onset time of the onset plots show the seconds of such times. The resample frequency must be a fixed duration (e.g. not a month).
- VDA plots: a legend of the symbols in the upper left corner of the plot, and a table of the results under it, with a row per fit: the release time at the Sun, the time its light reaches 1 AU (for comparison with observations at Earth), their error and the APL. Every fit line and its band start at inverse beta 0, where an open circle marks the release time at the Sun, also without reference times. The results were in the legends above the plots.
- The example events files (`examples/datetime_range_example.csv`, `datetime_range_bg_example.csv`, `datetime_range_only_example.csv` and `reference_times_example.csv`) have three Solar Orbiter events with a clear velocity dispersion, 2020-11-29, 2021-10-28 and 2022-01-20, whose VDA fits have no warnings with the default parameters. Their previous events (2021-05-22, 2021-11-09 and 2023-04-09) gave fits that are not physical. The example of the library uses the default viewings and grouped channels of the observer.
- The intermediate tables are not displayed by default (parameter `view_dfs` is `False`): they are shown with "Display the produced DataFrames" in the Views/Plots tab, or `tool.parameters.view_dfs = True`.
- A VDA fit of two onset points has no errors (`n/a`, NaN in the results) instead of errors of 0, and no band in the plot. A fit needs onsets at two different energies: onsets of one energy give no results instead of a fit of a vertical line.
- `vda.compute_vda()` gives results for every event (NaN without enough onset points, instead of no row for an event without onsets), and no longer keeps the results of the events of a previous run. The message of the events without enough onset points is printed with the results.
- The data step prints one line per event, e.g. `Event 1: data from 2021-10-28 14:00 to 2021-10-28 20:00, downloaded if they are not in the data folder`, instead of the messages and progress bars of the data loaders (e.g. "Creating dir", "Multi-index function only available…", "Unable to obtain … data!", which repeated the message of the tool for a spacecraft without data). If a loader fails, its messages are added to the error (`observers.quiet_loader`).
- The Solar Orbiter SPICE kernels are downloaded once, with one message and without progress bars, instead of again every 10 days. The tracebacks "Exception ignored ... signal only works in main thread" of their download in Jupyter, and the log of spiceypy when the tool is imported, are no longer shown.
- Labels: the parameters of the onset methods (e.g. "s: threshold in standard deviations above the background mean"), the time axis of the VDA plots ("Time (UTC)"), and the title of the onset plot (e.g. "HET protons"). The onset plot notes a zero background, whose levels cannot be shown in log scale.

### Fixed
- The links of the notebook text to the user guide go to its page on GitHub. They were broken in the folder downloaded with `spearhead-vda-notebooks`, which has no `docs` folder.
- The apparent path length (APL) is in AU: the path length $L = a \cdot c$ of the slope $a$ of the VDA fit. It was the path length divided by the distance of the observer from the Sun, too large by the factor 1 / distance (e.g. 1.76 instead of 1.42 AU for Solar Orbiter at 0.80 AU in the default event). The unit is shown in the printed results and the VDA plots.
- The energy ranges of the channels are read again in the Data step when the events or the spacecraft were changed after the parameters form was shown (e.g. an events file typed in the form). The energy ranges of the default event were used, which differ for SOHO EPHIN events before October 2017 and for Wind.
- Solar Orbiter without data on SOAR for an event (e.g. 2021-09-18) gives the error `No Solar Orbiter HET data from … to …`, as the other spacecraft, instead of `TypeError: list indices must be integers or slices, not str`. In the multi-spacecraft notebook, Solar Orbiter is then left out with this message.

### Upgrading from 1.1.0
- Update the installation: `pip install --upgrade "spearhead-vda[notebook]"` (from PyPI), or update the files of the tool for an installation from its folder. The only new dependency is `pooch`, which was already installed with seppy.
- The notebooks of the new version are downloaded with `spearhead-vda-notebooks --overwrite [folder]`; notebooks of 1.1.0 keep working.
- To see the intermediate tables as before, tick "Display the produced DataFrames" in the Views/Plots tab.
- The release times of resampled data are later by half the resample frequency (2.5 minutes for 5 minutes), and onsets close to their threshold can change (see Changed). The APL is unchanged where the onsets only move by half a bin.
- Data saved to a .pkl file by 1.1.0 or earlier with resampling are timestamped at the start of their bins (a message says so when they are loaded): download them again for the times at the centre of the bins.
- Otherwise the results are unchanged, except for SOHO EPHIN events before October 2017 and Wind events analysed with an events file or events changed after the form was shown (now with the energy ranges of their first event), and for fits of two points (errors `n/a` instead of 0).
- The APL of every analysis changes, except for an observer at 1 AU: multiply a previous APL by the distance of the observer from the Sun (AU) to compare. The release times are unchanged.
- Reference times are now entered as observed: a `(label, datetime)` reference time is observed at Earth, and is shown about 8 minutes earlier than before, at the Sun. Give the observer of the others, or their distance from the Sun.
- The results table has two new columns before `Release Time` (still the time at the spacecraft): `Release Time (Sun)` and `Release Time (1 AU)`. Code that reads the columns by name is unaffected. The printed results have the release time at the Sun first, and no "Extra Time" line (the light travel time is on the line of the time at the spacecraft).
- `analysis.fit_vda` returns the warnings of the fit in its second value, and no longer prints a message for two points.

## [1.1.0] - 2026-10-17

### Added
- Parameter `data_path`: the folder of the downloaded data files, by default `particle_data` in the current directory (as before). It is in the Data tab of the multi-event notebook and in the Event tab of the multi-spacecraft notebook, where it is the same for all the spacecraft. A relative folder is relative to the current directory when the data are downloaded, and `~` is the home folder. See the [user guide](docs/user_guide.md#parameters).

### Changed
- The downloaded Solar Orbiter data files are in the subfolder `solo` of the data folder, as those of the other spacecraft are in their own subfolders. Solar Orbiter data downloaded with previous versions (in `particle_data/l2`) are downloaded again once, and the old files can be deleted.

### Upgrading from 1.0.0
- Update the installation: `pip install --upgrade "spearhead-vda[notebook]"` (from PyPI), or update the files of the tool (download them again, or `git pull`) for an installation from its folder. The dependencies are unchanged.
- The notebooks of the new version are downloaded with `spearhead-vda-notebooks --overwrite [folder]`. Notebooks of 1.0.0 keep working, and show the data folder field too (only their text does not mention it).
- The data are still downloaded to `particle_data` in the current directory, except for Solar Orbiter, now in `particle_data/solo`.
- The results are unchanged.

## [1.0.0] - 2026-10-14

### Added
- The package `spearhead-vda` is published on [PyPI](https://pypi.org/project/spearhead-vda/): `pip install "spearhead-vda[notebook]"` installs the tool, and the command `spearhead-vda-notebooks [folder]` downloads its notebooks and example files (from the GitHub release of the installed version). The package metadata has classifiers and keywords, and the links of the README are absolute, for its page on PyPI.
- Notebook `vda_multi_spacecraft.ipynb`: the VDA of one event from several spacecraft, with the VDA fits of the spacecraft in one plot. As in the multi-event notebook, it has a parameters form (tabs Event, Spacecraft, Onsets and Views/Plots) and one cell per step; each spacecraft has its viewings, resample frequency and grouped channels, and the event, onset determination and selection, and plots are the same for all of them. A spacecraft without data for the event is left out. The VDA plot of each spacecraft can also be shown. Its interface is `spearhead.vda.notebook.VDA_multi_spacecraft_notebook`, and the combined plot is `views.plot_vda_combined`. See the [user guide](docs/user_guide.md#multiple-spacecraft).
- Notebook `vda_custom_points.ipynb`: the VDA fit of onset points given by the user, with an onset time, energy (single energy or channel range, MeV) and particle per point, in groups (e.g. one per spacecraft) with the distance of their observer from the Sun. The points are entered in its form (tabs Points and Views/Plots) or read from a points file (see `examples/custom_points_example.csv`), and the fits of the groups are plotted together, with an optional title. In scripts, `analysis.CustomPointsVDA` and `analysis.read_points_file`. See the [user guide](docs/user_guide.md#custom-points).
- Reference times in the VDA plot: optional (label, datetime) times, e.g. of the SXR peak or of type II and type III radio bursts, shown at inverse beta 0 with the fit line extended to it. In the notebook they are added in the Views/Plots tab ("Add Reference Time", or `tool.reference_times`), each shown in the plot of the event whose data range contains it; in scripts, `views.plot_vda(vda, event_no, references=[...])`.
- README: example VDA plots of the three notebooks (in `docs/images/`).
- Messages for the data servers that cannot be reached (e.g. a connection aborted while downloading): the parameters form and the data step say that the server may be unavailable, and to try again later.

### Changed
- The notebook `vda_tool.ipynb` is renamed `vda_multi_event.ipynb` (the VDA of one or more events from one spacecraft).
- The spacecraft of the data is named in the printed results and in the title of the VDA plot, e.g. `Event 1 (2021-10-28) | Solar Orbiter`.
- The inverse beta of a channel (`analysis.inverse_beta`, `analysis.mean_energy`) and the VDA fit (`analysis.fit_vda`) are functions shared by the VDA of the data and of custom points; the results are unchanged.

### Upgrading from 0.8.0
- Update the files of the tool (download them again, or `git pull`). The dependencies are unchanged: an installation from the folder of the tool (`pip install -e ".[notebook]"` or `pip install -r requirements.txt`) does not need to be updated. An installation from GitHub is updated with `pip install "spearhead-vda @ git+https://github.com/spearhead-he/VDA@v1.0.0"`.
- The tool can now also be installed from PyPI, in a new environment: `pip install "spearhead-vda[notebook]"`, and the notebooks downloaded with `spearhead-vda-notebooks [folder]` (see the [README](README.md#install-locally)).
- Open `vda_multi_event.ipynb` instead of `vda_tool.ipynb`. A copy of `vda_tool.ipynb` with your own changes keeps working.
- The results are unchanged.

## [0.8.0] - 2026-10-13

### Changed
- pyonset is installed from PyPI (`pyonset>=1.2.2`) instead of GitHub, so installing the tool from its folder (or from requirements.txt) no longer needs `git`.

### Upgrading from 0.7.0
- Update the installation once, from the folder of the tool: `pip install -r requirements.txt` (or `pip install -e ".[notebook]"`). This replaces pyonset from GitHub with pyonset 1.2.2 or newer from PyPI.
- The results are unchanged.

## [0.7.0] - 2026-10-13

### Changed
- VDA plot: the onsets of the electron channels are shown as stars, those of the proton channels as circles, with a legend entry for each. `vda_fits` records the particle of each point (`"particles"`).
- Tests: the offline tests are in one file per part of the tool (`tests/test_onsets.py`, `test_events.py`, `test_data.py`, `test_observers.py`, `test_onset_selection.py`, `test_vda_fit.py`, `test_parameters.py`, `test_notebook.py`) instead of `tests/test_units.py`, with new tests of the plots (`test_views.py`) and of the notebook interface (`test_notebook.py`: the parameters form, the background window and onset selection widgets, and the steps of the analysis), and the end-to-end test is `tests/test_end_to_end.py`. The end-to-end test, which downloads data, is marked `network`: `pytest -m "not network" tests/` runs only the offline tests. The commands to run the tests and update the baseline image are in [docs/development.md](docs/development.md).

### Removed
- Support for Python 3.10, which reaches its end of life in October 2026: the tool requires Python 3.11 or newer. The `notebook` extra requires ipykernel 7.4 or newer on all Python versions.

### Fixed
- The data of the sensors and viewings of an event are kept in time order also with a data gap at the start of one of them, as future pandas versions will not sort them by default (pandas deprecation warning).

### Upgrading from 0.6.0
- Python 3.11 or newer is required. With Python 3.10, create a new environment with a newer Python (e.g. `conda create --name vda python=3.12`) and install the tool in it, from the folder of the tool: `pip install -r requirements.txt`. With Python 3.11 or newer, the installation does not need to be updated.
- The results are unchanged. In the VDA plot, the onsets of the electron channels are stars.

## [0.6.0] - 2026-10-10

### Added
- Poisson-CUSUM onset determination method (`onset_method = "poisson_cusum"`), from pyonset, with the parameters `cusum_minutes` and `sigma_multiplier`. The onset plots show its background level, μd and the k and h parameters.
- Observers: the spacecraft of the analysis is the parameter `observer` (default `"solo"`, Solar Orbiter). The module `spearhead.vda.observers` describes each observer (sensors, particles, energy channels, viewings, default grouped channels), loads its data and gives its distance from the Sun.
- STEREO-A observer (`observer = "sta"`): IMPACT/HET (omni viewing) and IMPACT/SEPT (sun, asun, north and south viewings), loaded with seppy. Its distance from the Sun is from JPL Horizons. Its default grouped channels are the consecutive HET channels, in groups of three for protons and two for electrons (a single channel left at the end joins the previous group). The observer is selected in the Data tab of the parameters form; without data of the observer for the first event, the form shows an error in red and grays out the grouped channels.
- Parker Solar Probe observer (`observer = "psp"`): ISOIS/EPI-Hi HET, with the viewings A (sunward) and B (anti-sunward), from the 1-minute rates (`PSP_ISOIS-EPIHI_L2-HET-RATES60`) loaded with seppy. The proton channels are 3 to 11 (the others have no data), and the electron data are count rates divided by the channel width. Its distance from the Sun is from JPL Horizons.
- SOHO observer (`observer = "soho"`): ERNE-HED protons (`SOHO_ERNE-HED_L2-1MIN`, 13-64 MeV) and COSTEP-EPHIN electrons (level 2 data of the University of Kiel, E150 and E1300 with their energy ranges of the date; events on both sides of the change of the E1300 range, on 4 October 2017, raise an error), loaded with seppy, with the omni viewing. Its distance from the Sun is from JPL Horizons.
- Wind observer (`observer = "wind"`): the omnidirectional 3DP fluxes of the protons (SST Open, `WI_SOSP_3DP`) and electrons (SST Foil, `WI_SFSP_3DP`), loaded with seppy, with the energy widths of seppy (60% of the mean energy). Its distance from the Sun is from JPL Horizons.
- BepiColombo observer (`observer = "bepi"`): the SIXS-P level 3 cruise phase data of the SERPENTINE project, loaded with seppy, with the detector sides 0 to 2 as viewings (side 0 perpendicular to the Sun direction, sides 1 and 2 anti-sunward, from the SPICE attitude of the spacecraft). Its distance from the Sun is from JPL Horizons.
- Viewings of each sensor: each sensor uses the selected viewings it has (`VDA.sensor_viewings(sensor)`), also in the plots. A sensor of the grouped channels without any selected viewing raises an error.
- Saved data (`save_data_filepath`) include their observer, and loading the data of another observer raises an error. Data saved before are Solar Orbiter data.

### Changed
- The onset of the Sigma threshold method is the last point before the first of the `n` points above the threshold (previously the first point above it), as in the Poisson-CUSUM method of pyonset. The onsets, and so the release times, are one data point earlier (e.g. 5 minutes with the default resampling); the apparent path lengths are unchanged when all the onsets move by the same time.
- Setting `onset_method` sets `onset_method_parameters` to the defaults of the method, and the parameters form shows the parameters of the selected method. The onset methods are listed with their names (Sigma threshold, Poisson-CUSUM).
- Parameters of another onset method raise an error when the onsets are calculated.
- `viewings` and `channel_groups` default to `None`, replaced by the default viewings and grouped channels of the observer (unchanged for Solar Orbiter). Setting `observer` sets them to the defaults of the new observer.
- `seppy` (0.5.1 or newer) is a requirement again, for the data of the new spacecraft (it was already installed by pyonset).
- `VDA_parameters.AVAILABLE_SENSORS_PARTICLES`, `AVAILABLE_CHANNELS` and `AVAILABLE_VIEWINGS` are those of the observer.

### Fixed
- After the onset determination of a channel failed, the next channels of the event were determined without a background window.

### Removed
- The unused and not selectable Poisson-CUSUM bootstrap method (`_onset_detection_poisson_cusum_bootstrap`).
- The Solar Orbiter constants of `conf` (`AVAILABLE_SENSORS_PARTICLES`, `AVAILABLE_CHANNELS`, `AVAILABLE_VIEWINGS`, `DEFAULT_CHANNEL_GROUPS`) and of `VDA` (`RAW_FLUX_COLUMN`, `RAW_ENERGY_BINS_COLUMN`): they are attributes of `observers.SolarOrbiter`.

### Upgrading from 0.5.0
- Update the installation once, from the folder of the tool: `pip install -r requirements.txt` (or `pip install -e ".[notebook]"`). This installs seppy 0.5.1 or newer, used for the data of the new spacecraft.
- The Sigma threshold onsets, and so the release times, are one data point earlier than with 0.5.0 (e.g. 5 minutes with the default resampling); the apparent path lengths are unchanged when all the onsets move by the same time.
- The Solar Orbiter constants of `spearhead.vda.conf` (`AVAILABLE_SENSORS_PARTICLES`, `AVAILABLE_CHANNELS`, `AVAILABLE_VIEWINGS`, `DEFAULT_CHANNEL_GROUPS`) are attributes of `spearhead.vda.observers.SolarOrbiter`, or of the parameters (`VDA_parameters.AVAILABLE_SENSORS_PARTICLES` etc., for the selected observer).
- Data files saved with 0.5.0 (`save_data_filepath`) can still be loaded, with `observer = "solo"`: they are Solar Orbiter data, the only spacecraft of 0.5.0.

## [0.5.0] - 2026-10-09

### Added
- Installable Python package `spearhead-vda`, imported as `spearhead.vda`: `VDA`, `VDA_parameters` and `OnsetSelection` are imported from `spearhead.vda`, and the modules are `analysis` (previously `vda.py`), `conf` (previously `vda_tool_configuration.py`), `views` and `notebook` (previously `vda_notebook.py`). It is installed with `pip install -e ".[notebook]"` from the folder of the tool (or `pip install -r requirements.txt`), or `pip install "spearhead-vda @ git+https://github.com/spearhead-he/VDA@v0.5.0"`. The extras are `notebook` (Jupyter and ipywidgets) and `test`. The notebook also runs from the folder of the tool without installing the package. The version is `spearhead.vda.__version__`.
- `spearhead.vda.views` with the plots as functions that return matplotlib figures (`plot_bg`, `plot_onsets`, `plot_vda`), without displaying them. They are saved only when a filename is given.
- `VDA.select_onsets()`: the "Use all" onset selection, without widgets.
- `channel_groups` has the default grouped channels as its default value, so they are also used without the notebook. The channel groups widget shows the groups of `channel_groups`.

### Changed
- Wider background plot, so that its title and legend fit with large fonts.
- The widgets of `vda_views.py` (`VDA_nb_displayer`) are part of `VDA_notebook`. `VDA_notebook(vda)` can use an existing `VDA` object, e.g. for `background_selection()` and `onset_selection()` (the interactive onset selection) in another notebook.
- `bg_window`, `bg_window_source`, `bg_window_points`, `bg_window_warnings`, `check_bg_window`, `channel_energy_range`, `format_timedelta` and `vda_fits` of `VDA` are public (previously with a leading underscore).

### Fixed
- The notebook could stop responding (a cell never started) with ipykernel 7.0 to 7.3, mostly after the background window step. The `notebook` extra requires ipykernel 7.4 or newer on Python 3.11 and newer, and ipykernel 6 on Python 3.10.

### Removed
- `astrospice` and `seppy` from the requirements: they are not used by the tool (`seppy` is installed by `pyonset`).
- `VDA.plot()`, `VDA.plot_vda()` and `VDA.plot_bg_selection()`: replaced by `VDA.compute_vda()` and the functions of `spearhead.vda.views`.
- Parameter `default_channel_groups`: the default groups are the default value of `channel_groups`.
- `VDA_nb_displayer`: its widgets are part of `VDA_notebook`, its `construct_energies_df()` and `select_onsets()` are replaced by `VDA.construct_energies_df()` and `VDA.select_onsets()`.

### Upgrading from 0.4.0
- Install the package once, from the folder of the tool: `pip install -r requirements.txt` (or `pip install -e ".[notebook]"`). This also updates ipykernel in existing environments, which fixes the notebook hang.
- Imports: `from vda import VDA` becomes `from spearhead.vda import VDA`, `vda_tool_configuration` becomes `spearhead.vda.conf`, `vda_views` becomes `spearhead.vda.views` (now the plots) and `vda_notebook` becomes `spearhead.vda.notebook`.
- `vda.plot()` becomes `vda.compute_vda()` followed by `views.plot_vda(vda, event_no)` for each event (`views.plot_bg(vda, event_no)` instead of `vda.plot_bg_selection()`).
- `VDA_nb_displayer(vda).select_onsets()` becomes `vda.select_onsets()`. The background and interactive onset selection widgets are `VDA_notebook(vda).background_selection()` and `VDA_notebook(vda).onset_selection()`.
- `default_channel_groups` is replaced by `channel_groups`, whose default value has the same groups.

## [0.4.0] - 2026-10-06

### Added
- Notebook interface (`vda_notebook.py`): a form with all the parameters, in the tabs Events, Data, Energy channels, Onsets and Views/Plots, and one notebook cell per step of the analysis. The notebook has 6 code cells instead of 26.
- Energy range of each channel in the channel selection lists, instead of a separate energies table.
- The resample frequency is saved with the downloaded data, and used when the data are loaded.
- `viewings` is a list of viewing names, whose order is the priority of the "Use all" onset selection.
- Offline unit tests (`tests/test_units.py`), run before the end-to-end test. CI runs all the tests in `tests/`.

### Changed
- `VDA_parameters` is a dataclass: it can be printed, compared and copied with `dataclasses.replace`.
- Setting a removed or unknown parameter (e.g. a typo) raises an error, instead of being silently ignored. The errors of removed parameters explain their replacement.
- Loading and saving the data are set only with their file paths: an empty path means that the data are downloaded, or not saved.
- `OnsetSelection` enum for the onset selection method.
- The SPICE kernels information is no longer printed in the VDA step.
- `VDA.construct_energies_df` builds the energies table (it was in `vda_views.py`), and the widgets of `vda_views.py` are split into reusable parts.

### Removed
- "Custom List" onset selection, which was not implemented.
- Parameters `load_data`, `save_data` and `viewings_tt`.

### Upgrading from 0.3.0
- The notebook has a new layout: use the parameters form and the step cells. The analysis object is `tool.vda` and its parameters `tool.parameters`.
- `viewings_tt` is replaced by `viewings`, e.g. `["sun", "north"]`.
- `load_data` / `save_data` are replaced by `load_data_filepath` / `save_data_filepath` alone.
- `onset_selection` accepts `OnsetSelection.USE_ALL` (0) and `OnsetSelection.INTERACTIVE` (1).

## [0.3.0] - 2026-09-30

### Added
- Events file whose type is deduced from its columns: datetime ranges (`Start Time, End Time`, optionally with `BG End` or `BG Start, BG End`) or reference times (`Reference Time`). Example files for every layout in `examples/`.
- Several events from datetime ranges entered in the notebook, without an events file.
- Background window of each event in the events table: from the events file, or a default in minutes after the start of the data. It can be checked and changed per event interactively, and saved with `save_times()` for later runs.
- Warnings for background windows outside the data or with less than 3 points.
- Interactive onset selection one channel at a time, with event and channel dropdowns and previous / next buttons, starting from the "Use all" selection.
- `compute_vda()` and `print_results()`: the VDA results of every event, also without plots.
- New energy channel grouping selection.
- Tests with pytest and coverage, citation file, Python 3.14 support.

### Changed
- The background window is set with times instead of point indices, so it no longer depends on the resample frequency.
- Clearer background, onset and VDA plots.
- Code restructured to remove duplication.

### Fixed
- The VDA fit depended on the local timezone of the computer: onsets around a daylight saving time change gave wrong release times and path lengths.
- "Use all" used the last viewing with an onset instead of the first one.
- Reference times input, plotting of all the events, and the order of the particles.
- The background window of the datetime ranges file was ignored if the onset method parameters cell was not run.

### Removed
- Parameters `input_type`, `date_range_filepath`, `reference_times_filepath`, `date_start`, `date_end`, and the `bg_start` / `bg_end` onset method parameters.

### Upgrading from 0.2.0
- The input type selection is replaced by a single events file (`input_filepath`), whose type is deduced from its columns. Without a file, `date_ranges` (list of (start, end) datetimes) replaces `date_start` / `date_end`.
- The background window is given in the events file or as a default in minutes after the start of the data, instead of point indices (`bg_start` / `bg_end` of the onset method).
- Code that sets a removed parameter gets an error explaining its replacement.

## [0.2.0] - 2025-06-26

### Added
- `VDA`, `VDA_parameters` and `VDA_nb_displayer` classes, used by the notebook cells.
- Interactive onset selection.
- Option to not save the VDA plot.

### Fixed
- Limits of the background plot.

## [0.1.0-alpha] - 2024-12-13

First release: notebook for the Velocity Dispersion Analysis of Solar Energetic Particle events with Solar Orbiter EPD (HET, EPT) data.

[Unreleased]: https://github.com/spearhead-he/VDA/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/spearhead-he/VDA/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/spearhead-he/VDA/compare/v0.8.0...v1.0.0
[0.8.0]: https://github.com/spearhead-he/VDA/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/spearhead-he/VDA/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/spearhead-he/VDA/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/spearhead-he/VDA/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/spearhead-he/VDA/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/spearhead-he/VDA/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/spearhead-he/VDA/compare/v0.1.0-alpha...v0.2.0
[0.1.0-alpha]: https://github.com/spearhead-he/VDA/releases/tag/v0.1.0-alpha
