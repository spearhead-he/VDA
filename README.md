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
  - [What's new in v0.5.0](#whats-new-in-v050)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
  - [Acknowledgement](#acknowledgement)

## About

The VDA tool helps in the automation of Velocity Dispersion Analysis (VDA) of one or multiple Solar Energetic Particle (SEP) events. The events are provided to the tool as datetime ranges or reference times, from a file or entered in the Notebook. The user can parameterize the given Notebook and control the spacecraft whose data are used (Solar Orbiter, STEREO-A, Parker Solar Probe, SOHO, Wind or BepiColombo), which particle species are used, the sensor from which they are detected, the viewings to be considered, and the background window of each event.

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
6. Install the tool and its dependencies with `pip install -e ".[notebook]"` (or, equivalently, `pip install -r requirements.txt`)
7. Open the Jupyter Notebook by running `jupyter-lab vda_tool.ipynb`

## How to use

Open the Notebook `vda_tool.ipynb` and run its cells in order: a setup cell, a form with all the parameters (tabs Events, Data, Energy channels, Onsets and Views/Plots), and one cell per step of the analysis:

```python
from spearhead.vda.notebook import VDA_notebook
tool = VDA_notebook()

tool.parameters_form()       # all the parameters
tool.run_data()              # events, data download and grouping of the energy channels
tool.background_selection()  # background window of each event
tool.run_onsets()            # onsets, and their selection
tool.run_vda()               # results and plots
```

The events are given with a .csv file (datetime ranges or reference times, see the [examples](examples/)) or entered in the form. See the [user guide](docs/user_guide.md) for details.

## What's new in v0.5.0

- The tool is an installable Python package, `spearhead.vda`, which can also be used in scripts and other notebooks (see [Using the VDA tool without the notebook](docs/library.md)).
- The plots are returned as matplotlib figures (`spearhead.vda.views`), saved only when a filename is given.
- The background window and interactive onset selection widgets can be used in other notebooks with `VDA_notebook(vda)`.
- The default grouped channels are also used without the notebook.
- Fixed: the Notebook could stop responding with ipykernel 7.0 to 7.3.

See the [changelog](CHANGELOG.md) for all the changes, and how to upgrade from previous versions.

## Documentation

- [User guide](docs/user_guide.md): the Notebook, the spacecraft and their instruments, the parameters, the events files, background windows, onset selection and results
- [Using the VDA tool without the notebook](docs/library.md): the analysis code in scripts
- [Changelog](CHANGELOG.md)

## Contributing

Contributions to this tool are very much welcome and encouraged! Contributions can take the form of [issues](https://github.com/spearhead-he/VDA/issues) to report bugs and request new features or [pull requests](https://github.com/spearhead-he/VDA/pulls) to submit new code. 

If you don't have a GitHub account, you can [sign-up for free here](https://github.com/signup), or you can also reach out to us with feedback by sending an email to jan.gieseler@utu.fi.

## Acknowledgement

<img align="right" height="80px" src="https://github.com/user-attachments/assets/28c60e00-85b4-4cf3-a422-6f0524c42234"> 
<a href="https://spearhead-he.eu"><img align="right" height="80px" src="https://github.com/user-attachments/assets/854d45ef-8b25-4a7b-9521-bf8bc364246e"></a>

This tool is developed within the [SPEARHEAD (*SPEcification, Analysis & Re-calibration of High Energy pArticle Data*)](https://spearhead-he.eu/) project. SPEARHEAD has received funding from the European Union’s Horizon Europe programme under grant agreement No 101135044. 

The tool reflects only the authors’ view and the European Commission is not responsible for any use that may be made of the information it contains.
