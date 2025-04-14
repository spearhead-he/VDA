[![DOI](https://zenodo.org/badge/901860658.svg)](https://doi.org/10.5281/zenodo.14441053)

# SPEARHEAD VDA tool

- [SPEARHEAD VDA tool](#spearhead-vda-tool)
  - [About](#about)
  - [How to install](#how-to-install)
  - [How to use](#how-to-use)
  - [Contributing](#contributing)
  - [Acknowledgement](#acknowledgement)

## About

The VDA tool helps in the automation of Velocity Dispersion Analysis (VDA) of one or multiple Solar Energetic Particle (SEP) events. Each event is provided to the tool as a point in time. The user can parameterize the given Notebook and control which particle species are used, the sensor from which they are detected, and the viewings to be considered.

The tool utilizes the Pandas module and generates multiple DataFrames during its execution. The final output of the tool is a plot of the VDA analysis for each inputted event.

*Tool is still under active development and its results should be handled with caution.* 

*Tested in Ubuntu 22.04 with Python version 3.10.12, and MacOS 15.1.1 with Python 3.10.16 and 3.12.8*

## How to install

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
- User Inputs
- Run

The user should run the cell(s) of the first section and then follow the instructions iside the Notebook to properly fill the input forms. The cells of the "Run" section can then be run without changing anything.

## Contributing

Contributions to this tool are very much welcome and encouraged! Contributions can take the form of [issues](https://github.com/spearhead-he/VDA/issues) to report bugs and request new features or [pull requests](https://github.com/spearhead-he/VDA/pulls) to submit new code. 

If you don't have a GitHub account, you can [sign-up for free here](https://github.com/signup), or you can also reach out to us with feedback by sending an email to jan.gieseler@utu.fi.

## Acknowledgement

<img align="right" height="80px" src="https://github.com/user-attachments/assets/28c60e00-85b4-4cf3-a422-6f0524c42234"> 
<img align="right" height="80px" src="https://github.com/user-attachments/assets/854d45ef-8b25-4a7b-9521-bf8bc364246e"> 

This tool is developed within the SPEARHEAD (*SPEcification, Analysis & Re-calibration of High Energy pArticle Data*) project. SPEARHEAD has received funding from the European Union’s Horizon Europe programme under grant agreement No 101135044. 

The tool reflects only the authors’ view and the European Commission is not responsible for any use that may be made of the information it contains.
