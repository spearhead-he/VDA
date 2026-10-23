[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.14441053-blue)](https://doi.org/10.5281/zenodo.14441053)
[![Python versions](https://img.shields.io/badge/python-3.11_--_3.14-blue)]()
[![pytest](https://github.com/spearhead-he/VDA/actions/workflows/pytest.yml/badge.svg?branch=main)](https://github.com/spearhead-he/VDA/actions/workflows/pytest.yml)
[![codecov](https://codecov.io/github/spearhead-he/VDA/graph/badge.svg?token=GH3JBH0EGW)](https://codecov.io/github/spearhead-he/VDA)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![website](https://img.shields.io/badge/Project%20Website-blue)](https://spearhead-he.eu)

# SPEARHEAD VDA tool

- [SPEARHEAD VDA tool](#spearhead-vda-tool)
  - [About](#about)
  - [How to install](#how-to-install)
  - [How to use](#how-to-use)
  - [Example output](#example-output)
  - [What's new in v1.2.0](#whats-new-in-v120)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
  - [Acknowledgement](#acknowledgement)

## About

The VDA tool helps in the automation of Velocity Dispersion Analysis (VDA) of one or multiple Solar Energetic Particle (SEP) events. The events are provided to the tool as datetime ranges or reference times, from a file or entered in the Notebook. The user can parameterize the given Notebook and control the spacecraft whose data are used (Solar Orbiter, STEREO-A, Parker Solar Probe, SOHO, Wind or BepiColombo), which particle species are used, the sensor from which they are detected, the viewings to be considered, and the background window of each event.

The tool has three notebooks: `vda_multi_event.ipynb`, the VDA of one or more events from one spacecraft, `vda_multi_spacecraft.ipynb`, the VDA of one event from several spacecraft, with their VDA fits in one plot, and `vda_custom_points.ipynb`, the VDA fit of onset points given by the user, in groups such as one per spacecraft.

The tool utilizes the Pandas module and generates multiple DataFrames during its execution. The final output of the tool is the release time and apparent path length of each inputted event, with a plot of its VDA analysis.

*Tested in Ubuntu with Python 3.11 to 3.14, and MacOS 15.1.1 with Python 3.12.8*

## How to install

### Access online (JupyterHub)

You can access the complete VDA tool online on the JupyterHub server of the SPEARHEAD project. All you need is a (free) [GitHub account](https://github.com/signup) for verification. [Access the tool by opening this link!](https://jupyterhub.spearhead-he.eu/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Fspearhead-he%2FVDA&urlpath=lab%2Ftree%2FVDA%2Fvda_multi_event.ipynb&branch=main)

### Install locally

1. This tool requires a recent Python (>=3.11) installation. [Following SunPy's approach, we recommend installing Python via miniforge (click for instructions).](https://docs.sunpy.org/en/stable/tutorial/installation.html#installing-python)
2. Open a terminal or the miniforge prompt, create a new virtual environment (e.g., `conda create --name vda python=3.12`) and activate it (e.g., `conda activate vda`).
3. Install the tool and its dependencies from [PyPI](https://pypi.org/project/spearhead-vda/) with `pip install "spearhead-vda[notebook]"`
4. Download the notebooks and their example files into a folder of your choice (here `vda`) with `spearhead-vda-notebooks vda`, and move to it with `cd vda`
5. Open a Jupyter Notebook by running `jupyter-lab vda_multi_event.ipynb` (or `jupyter-lab vda_multi_spacecraft.ipynb`, `jupyter-lab vda_custom_points.ipynb`)

### Install from the repository

To use the latest code of the repository or change it:

1. [Download this file](https://github.com/spearhead-he/VDA/archive/refs/heads/main.zip) and extract to a folder of your choice (or clone the repository [https://github.com/spearhead-he/VDA](https://github.com/spearhead-he/VDA) if you know how to use `git`).
2. Open a terminal or the miniforge prompt and move to the directory where the code is.
3. Create a new virtual environment (e.g., `conda create --name vda python=3.12`) and activate it (e.g., `conda activate vda`).
4. Install the tool and its dependencies with `pip install -e ".[notebook]"` (or, equivalently, `pip install -r requirements.txt`)
5. Open the Jupyter Notebook by running `jupyter-lab vda_multi_event.ipynb` (or `jupyter-lab vda_multi_spacecraft.ipynb`, `jupyter-lab vda_custom_points.ipynb`)

## How to use

Open the Notebook `vda_multi_event.ipynb` and run its cells in order: a setup cell, a form with all the parameters (tabs Events, Data, Energy channels, Onsets and Views/Plots), and one cell per step of the analysis:

```python
from spearhead.vda.notebook import VDA_notebook
tool = VDA_notebook()

tool.parameters_form()       # all the parameters
tool.run_data()              # events, data download and grouping of the energy channels
tool.background_selection()  # background window of each event
tool.run_onsets()            # onsets, and their selection
tool.run_vda()               # results and plots
```

The events are given with a .csv file (datetime ranges or reference times, see the [examples](https://github.com/spearhead-he/VDA/tree/main/examples/)) or entered in the form. See the [user guide](https://github.com/spearhead-he/VDA/blob/main/docs/user_guide.md) for details.

The Notebook `vda_multi_spacecraft.ipynb` has the same cells, with `VDA_multi_spacecraft_notebook` instead of `VDA_notebook`, for one event and several spacecraft (see the [user guide](https://github.com/spearhead-he/VDA/blob/main/docs/user_guide.md#multiple-spacecraft)). The Notebook `vda_custom_points.ipynb` has a form with the onset points (or a points file, see the [example](https://github.com/spearhead-he/VDA/blob/main/examples/custom_points_example.csv)) and the VDA step (see the [user guide](https://github.com/spearhead-he/VDA/blob/main/docs/user_guide.md#custom-points)).

## Example output

The VDA plots of the notebooks with their default parameters and example data.

`vda_multi_event.ipynb`: the default event (2021-10-28) from Solar Orbiter, with the onsets of the proton channels as circles and of the electron channels as stars, and the release time and apparent path length (APL) of the VDA fit.

<img src="https://raw.githubusercontent.com/spearhead-he/VDA/main/docs/images/vda_multi_event.png" alt="VDA plot of the 2021-10-28 event from Solar Orbiter" width="700">

`vda_multi_spacecraft.ipynb`: the same event from Solar Orbiter, STEREO-A and SOHO, with the VDA fit of each spacecraft extended to inverse beta 0.

<img src="https://raw.githubusercontent.com/spearhead-he/VDA/main/docs/images/vda_multi_spacecraft.png" alt="VDA fits of the 2021-10-28 event from Solar Orbiter, STEREO-A and SOHO" width="700">

`vda_custom_points.ipynb`: the example points of the 2024-11-21 event (`examples/custom_points_example.csv`) seen by SOHO, STEREO-A and Solar Orbiter, with reference times at inverse beta 0.

<img src="https://raw.githubusercontent.com/spearhead-he/VDA/main/docs/images/vda_custom_points.png" alt="VDA fits of the custom points of the 2024-11-21 event, with reference times" width="700">

## What's new in v1.2.0

- The release time at the Sun and the time its light reaches 1 AU are given in the printed results, the results table and the VDA plots. The VDA plots have a legend of the symbols and a table of the results, and every fit starts at inverse beta 0.
- The apparent path length (APL) is in AU. It was too large by the factor 1 / distance of the observer from the Sun.
- Resampled data are timestamped at the centre of their bins, so the release times are half a bin later (2.5 minutes for the default 5 minutes).
- Each reference time has an observer (Earth, a spacecraft, or its distance from the Sun in AU), and is shown at the Sun.
- Warnings for VDA fits that are not physical: an APL that is not positive or is shorter than the distance of the observer from the Sun, or a release after the earliest onset.
- `vda.run()` runs the analysis in scripts, from the events to the VDA fits.
- Clearer and quieter messages: errors before the download, one line per event in the data step, and a clear error for Solar Orbiter without data on SOAR.

See the [changelog](https://github.com/spearhead-he/VDA/blob/main/CHANGELOG.md) for all the changes, and how to upgrade from previous versions.

## Documentation

- [User guide](https://github.com/spearhead-he/VDA/blob/main/docs/user_guide.md): the Notebooks, the spacecraft and their instruments, the parameters, the events files, background windows, onset selection, results, multiple spacecraft and custom points
- [Using the VDA tool without the notebook](https://github.com/spearhead-he/VDA/blob/main/docs/library.md): the analysis code in scripts
- [Development](https://github.com/spearhead-he/VDA/blob/main/docs/development.md): installation for development, the tests and the baseline image of the VDA plot
- [Changelog](https://github.com/spearhead-he/VDA/blob/main/CHANGELOG.md)

## Contributing

Contributions to this tool are very much welcome and encouraged! Contributions can take the form of [issues](https://github.com/spearhead-he/VDA/issues) to report bugs and request new features or [pull requests](https://github.com/spearhead-he/VDA/pulls) to submit new code. 

If you don't have a GitHub account, you can [sign-up for free here](https://github.com/signup), or you can also reach out to us with feedback by sending an email to jan.gieseler@utu.fi.

## Acknowledgement

<img align="right" height="80px" src="https://github.com/user-attachments/assets/28c60e00-85b4-4cf3-a422-6f0524c42234"> 
<a href="https://spearhead-he.eu"><img align="right" height="80px" src="https://github.com/user-attachments/assets/854d45ef-8b25-4a7b-9521-bf8bc364246e"></a>

This tool is developed within the [SPEARHEAD (*SPEcification, Analysis & Re-calibration of High Energy pArticle Data*)](https://spearhead-he.eu/) project. SPEARHEAD has received funding from the European Union’s Horizon Europe programme under grant agreement No 101135044. 

The tool reflects only the authors’ view and the European Commission is not responsible for any use that may be made of the information it contains.
