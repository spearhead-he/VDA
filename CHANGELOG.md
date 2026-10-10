# Changelog

All notable changes of the SPEARHEAD VDA tool. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- Poisson-CUSUM onset determination method (`onset_method = "poisson_cusum"`), from pyonset, with the parameters `cusum_minutes` and `sigma_multiplier`. The onset plots show its background level, μd and the k and h parameters.

- Observers: the spacecraft of the analysis is the parameter `observer` (for now only `"solo"`, Solar Orbiter). The module `spearhead.vda.observers` describes each observer (sensors, particles, energy channels, viewings, default grouped channels), loads its data and gives its distance from the Sun.
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
- `seppy` (0.5.1 or newer) is a requirement again, for the STEREO-A data (it was already installed by pyonset).
- `VDA_parameters.AVAILABLE_SENSORS_PARTICLES`, `AVAILABLE_CHANNELS` and `AVAILABLE_VIEWINGS` are those of the observer.

### Fixed
- After the onset determination of a channel failed, the next channels of the event were determined without a background window.

### Removed
- The unused and not selectable Poisson-CUSUM bootstrap method (`_onset_detection_poisson_cusum_bootstrap`).
- The Solar Orbiter constants of `conf` (`AVAILABLE_SENSORS_PARTICLES`, `AVAILABLE_CHANNELS`, `AVAILABLE_VIEWINGS`, `DEFAULT_CHANNEL_GROUPS`) and of `VDA` (`RAW_FLUX_COLUMN`, `RAW_ENERGY_BINS_COLUMN`): they are attributes of `observers.SolarOrbiter`.

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

[Unreleased]: https://github.com/spearhead-he/VDA/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/spearhead-he/VDA/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/spearhead-he/VDA/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/spearhead-he/VDA/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/spearhead-he/VDA/compare/v0.1.0-alpha...v0.2.0
[0.1.0-alpha]: https://github.com/spearhead-he/VDA/releases/tag/v0.1.0-alpha
