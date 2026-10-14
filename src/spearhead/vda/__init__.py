"""Velocity Dispersion Analysis (VDA) of Solar Energetic Particle events.

The data are from one spacecraft (observer) at a time: Solar Orbiter, STEREO-A, Parker Solar Probe, SOHO, Wind or
BepiColombo.

Modules:
- analysis: the VDA class, which performs the analysis, and CustomPointsVDA, the VDA fit of onset points given by
  the user
- conf: VDA_parameters, the parameters of the analysis, and OnsetSelection
- observers: the spacecraft, their instruments and data loading
- views: the plots, as matplotlib figures
- notebook: the interfaces of the notebooks with their widgets, VDA_notebook (multiple events),
  VDA_multi_spacecraft_notebook and VDA_custom_points_notebook (needs the notebook dependencies)
- notebooks: download of the notebooks and their example files (the command spearhead-vda-notebooks)
"""
from .analysis import VDA
from .conf import OnsetSelection, VDA_parameters

__version__ = "1.0.0"

__all__ = ["VDA", "VDA_parameters", "OnsetSelection", "__version__"]
