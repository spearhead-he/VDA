# Changelog

All notable changes of the SPEARHEAD VDA tool. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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

[0.4.0]: https://github.com/spearhead-he/VDA/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/spearhead-he/VDA/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/spearhead-he/VDA/compare/v0.1.0-alpha...v0.2.0
[0.1.0-alpha]: https://github.com/spearhead-he/VDA/releases/tag/v0.1.0-alpha
