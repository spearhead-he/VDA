"""Velocity Dispersion Analysis (VDA) of Solar Energetic Particle events with Solar Orbiter EPD data.

Modules:
- analysis: the VDA class, which performs the analysis
- conf: VDA_parameters, the parameters of the analysis, and OnsetSelection
- views: the plots, as matplotlib figures
- notebook: VDA_notebook, the interface of the notebook with its widgets (needs the notebook dependencies)
"""
from .analysis import VDA
from .conf import OnsetSelection, VDA_parameters

__all__ = ["VDA", "VDA_parameters", "OnsetSelection"]
