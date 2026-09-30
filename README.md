[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.14441053-blue)](https://doi.org/10.5281/zenodo.14441053)
[![Python versions](https://img.shields.io/badge/python-3.10_--_3.14-blue)]()
[![pytest](https://github.com/spearhead-he/VDA/actions/workflows/pytest.yml/badge.svg?branch=main)](https://github.com/spearhead-he/VDA/actions/workflows/pytest.yml)
[![codecov](https://codecov.io/github/spearhead-he/VDA/graph/badge.svg?token=GH3JBH0EGW)](https://codecov.io/github/spearhead-he/VDA)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![website](https://img.shields.io/badge/Project%20Website-blue)](https://spearhead-he.eu)

# SPEARHEAD VDA tool

- [SPEARHEAD VDA tool](#spearhead-vda-tool)
  - [About](#about)
  - [How to install](#how-to-install)
  - [How to use](#how-to-use)
    - [Events](#events)
    - [Background window](#background-window)
    - [Onset selection](#onset-selection)
    - [Results](#results)
    - [Changes from v0.2.0](#changes-from-v020)
  - [Contributing](#contributing)
  - [Acknowledgement](#acknowledgement)

## About

The VDA tool helps in the automation of Velocity Dispersion Analysis (VDA) of one or multiple Solar Energetic Particle (SEP) events. The events are provided to the tool as datetime ranges or reference times, from a file or entered in the Notebook. The user can parameterize the given Notebook and control which particle species are used, the sensor from which they are detected, the viewings to be considered, and the background window of each event.

The tool utilizes the Pandas module and generates multiple DataFrames during its execution. The final output of the tool is the release time and apparent path length of each inputted event, with a plot of its VDA analysis.

*Tool is still under active development and its results should be handled with caution.* 

*Tested in Ubuntu 22.04 with Python version 3.10.12, and MacOS 15.1.1 with Python 3.10.16 and 3.12.8*

## How to install

### Access online (JupyterHub)

You can access the complete VDA tool online on the JupyterHub server of the SPEARHEAD project. All you need is a (free) [GitHub account](https://github.com/signup) for verification. [Access the tool by opening this link!](https://jupyterhub.spearhead-he.eu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Fspearhead-he%2FVDA&urlpath=lab%2Ftree%2FVDA%2Fvda_tool.ipynb&branch=main)

### Install locally

1. This tool requires a recent Python (>=3.10) installation. [Following SunPy's approach, we recommend installing Python via miniforge (click for instructions).](https://docs.sunpy.org/en/stable/tutorial/installation.html#installing-python)
2. [Download this file](https://github.com/spearhead-he/VDA/archive/refs/heads/main.zip) and extract to a folder of your choice (or clone the repository [https://github.com/spearhead-he/VDA](https://github.com/spearhead-he/VDA) if you know how to use `git`).
3. Open a terminal or the miniforge prompt and move to the directory where the code is.
4. Create a new virtual environment (e.g., `conda create --name vda python=3.12`) and activate it (e.g., `conda activate vda`).
5. If you **don't** have `git` installed (try executing it), install it with `conda install conda-forge::git`.
6. Install the Python dependencies from the *requirements.txt* file with `pip install -r requirements.txt`
7. Open the Jupyter Notebook by running `jupyter-lab vda_tool.ipynb`

## How to use

The Notebook is separated into three main sections:
- Imports & Setup
- Parameterize
- VDA

The user should run the cell(s) of the first section and then follow the instructions inside the Notebook to properly fill the input forms. The cells of the "VDA" section can then be run without changing anything.

### Events

The events are given with a .csv file, or as datetime ranges entered in the Notebook ("Add Event" button) if no file is given. The first column of the file is the event number, and the type of the file is deduced from the names of the other columns:

| Columns after the event number | Meaning | Example |
|---|---|---|
| `Start Time, End Time` | Datetime range of the data. The default background window is used | [examples/datetime_range_only_example.csv](examples/datetime_range_only_example.csv) |
| `Start Time, BG End, End Time` | As above, with the background window from `Start Time` to `BG End` | [examples/datetime_range_example.csv](examples/datetime_range_example.csv) |
| `Start Time, BG Start, BG End, End Time` | As above, with the background window from `BG Start` to `BG End`. Empty cells use the default background window | [examples/datetime_range_bg_example.csv](examples/datetime_range_bg_example.csv) |
| `Reference Time` | Reference datetime of the event. The datetime range of the data is set with the hours prior to and after it | [examples/reference_times_example.csv](examples/reference_times_example.csv) |

The datetimes can be in any format supported by `pandas.to_datetime` (e.g. `2024-12-31 00:00:00`).

### Background window

Each event has a background window, used for the onset determination:
- the background window of the events file, if given,
- otherwise the default background window, in minutes after the start of the event's data (slider in the Notebook).

The background window of each event can be checked and changed in the Notebook with the event dropdown and the "Background" slider, or set back to the default one. The chosen windows can be saved with `vda.save_times("path.csv")`, and the saved file can be used as the events file of later runs.

### Onset selection

- Use all: for each grouped energy channel, the onset of the first viewing (in the order of the viewings) with a determined onset is used.
- Interactive: the determined onsets are shown one channel at a time, chosen with the event and channel dropdowns or the "Previous" / "Next" buttons, and the viewing whose onset is used is selected per channel (or none, to leave the channel out). The selection starts from the viewings of "Use all".

### Results

For each event, the release time, the extra time (light travel time from the Sun to the spacecraft) and the apparent path length (APL) are printed, followed by the VDA plot, which is also saved as a .png file. The results of all the events are stored in the `vda.results` table. For many events, `vda.compute_vda()` followed by `vda.print_results()` gives the results without creating the plots.

### Changes from v0.2.0

- The input type selection is replaced by a single events file (`input_filepath`), whose type is deduced from its columns. Without a file, `date_ranges` (list of (start, end) datetimes) replaces `date_start` / `date_end`.
- The background window is given in the events file or as a default in minutes after the start of the data, instead of point indices (`bg_start` / `bg_end` of the onset method).
- Code that sets a removed parameter gets an error explaining its replacement.

## Contributing

Contributions to this tool are very much welcome and encouraged! Contributions can take the form of [issues](https://github.com/spearhead-he/VDA/issues) to report bugs and request new features or [pull requests](https://github.com/spearhead-he/VDA/pulls) to submit new code. 

If you don't have a GitHub account, you can [sign-up for free here](https://github.com/signup), or you can also reach out to us with feedback by sending an email to jan.gieseler@utu.fi.

## Acknowledgement

<img align="right" height="80px" src="https://github.com/user-attachments/assets/28c60e00-85b4-4cf3-a422-6f0524c42234"> 
<a href="https://spearhead-he.eu"><img align="right" height="80px" src="https://github.com/user-attachments/assets/854d45ef-8b25-4a7b-9521-bf8bc364246e"></a>

This tool is developed within the [SPEARHEAD (*SPEcification, Analysis & Re-calibration of High Energy pArticle Data*)](https://spearhead-he.eu/) project. SPEARHEAD has received funding from the European Union’s Horizon Europe programme under grant agreement No 101135044. 

The tool reflects only the authors’ view and the European Commission is not responsible for any use that may be made of the information it contains.
