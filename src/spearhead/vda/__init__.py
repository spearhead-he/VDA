"""Velocity Dispersion Analysis (VDA) of Solar Energetic Particle events.

The data are from one spacecraft (observer) at a time: Solar Orbiter, STEREO-A, Parker Solar Probe, SOHO, Wind or
BepiColombo.

Modules:
- analysis: the VDA class, which performs the analysis
- conf: VDA_parameters, the parameters of the analysis, and OnsetSelection
- observers: the spacecraft, their instruments and data loading
- views: the plots, as matplotlib figures
- notebook: VDA_notebook, the interface of the notebook with its widgets (needs the notebook dependencies)
"""
from .analysis import VDA
from .conf import OnsetSelection, VDA_parameters

__version__ = "0.6.0"

__all__ = ["VDA", "VDA_parameters", "OnsetSelection", "__version__"]
