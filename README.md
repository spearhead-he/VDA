# SPEARHEAD VDA tool

This tool is developed for the SPEARHEAD (*SPEcification, Analysis & Re-calibration of High Energy pArticle Data*) project. SPEARHEAD has received funding from the European Union’s Horizon Europe programme under grant agreement No 101135044.

**Table of Contents**
- [SPEARHEAD VDA tool](#spearhead-vda-tool)
  - [About the tool](#about-the-tool)
  - [How to install](#how-to-install)
  - [How to use](#how-to-use)

## About the tool

The VDA tool helps in the automation of the VDA analysis of one or multiple events. Each event is provided to the tool as a point in time. The user can parameterize the given notebook and control what particle species are used, the sensor from which they are detected, and the viewings to be considered.

The tool utilizes the pandas module and generates multiple DataFrames during its execution. The final output of the tool is a plot of the VDA analysis for each inputted event.

*Tool is tested in Ubuntu 22.04 with python version 3.10.12*

## How to install

1. Download/Clone the contents of this repository.
2. Move to the directory where the code is (unzip the downloaded file, if needed).
3. *(Optional)* Create a python virtual environment (e.g. `python -m venv venv_vda_tool`) and activate it (e.g. `source venv_vda_tool/bin/activate`)
4. Install the python dependencies from the *requirements.txt* file `pip install -r requirements.txt`

## How to use

The notebook is separated into 3 main sections:
- Imports & Setup
- User Inputs
- Run

The user should run the cell(s) of the frst section and then follow the instructions iside the notebook to properly fill the input forms. The cells of the "Run" section can then be run automatically.