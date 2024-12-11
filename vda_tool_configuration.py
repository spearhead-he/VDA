from ipywidgets import widgets
from os import getcwd

############### Widgets ###############
WIDGETS_LAYOUT = widgets.Layout(width="auto")
WIDGETS_STYLE = {"description_width": "initial"}

############### Reference Times DF ###############
REF_TIME_COLNAME = "Reference Time"
START_TIME_COLNAME = "Start Time"
END_TIME_COLNAME = "End Time"

############### Particle Data ###############
DATA_PATH = f"{getcwd()}/particle_data"
PROTON_COLUMN_PREFIX = "H_Flux"
ELECTRON_COLUMN_PREFIX = "Electron_Flux"

AVAILABLE_SENSORS = ["het", "ept"]
AVAILABLE_PARTICLES = ["protons", "electrons"]
AVAILABLE_VIEWINGS = ["sun", "asun", "north", "south", "omni"]

############### Particle Data ###############
AVAILABLE_ONSET_METHODS = {
    "sigma": {
        "s": {
            "widget": widgets.IntSlider,
            "widget_params": {
                "value": 3,
                "min": 1,
                "max": 5,
                "step": 1,
                "description": "s:",
                "disabled": False,
                "style": WIDGETS_STYLE,
                "layout": WIDGETS_LAYOUT
            }
        },
        "n": {
            "widget": widgets.IntSlider,
            "widget_params": {
                "value": 3,
                "min": 1,
                "max": 5,
                "step": 1,
                "description": "n:",
                "disabled": False,
                "style": WIDGETS_STYLE,
                "layout": WIDGETS_LAYOUT
            }
        },
        "bg_start": {
            "widget": widgets.IntSlider,
            "widget_params": {
                "value": 0,
                "min": 0,
                "max": 100,
                "step": 1,
                "description": "Point index to start the background sampling:",
                "disabled": False,
                "style": WIDGETS_STYLE,
                "layout": WIDGETS_LAYOUT
            }
        },
        "bg_end": {
            "widget": widgets.IntSlider,
            "widget_params": {
                "value": 12,
                "min": 0,
                "max": 100,
                "step": 1,
                "description": "Point index to end the background sampling:",
                "disabled": False,
                "style": WIDGETS_STYLE,
                "layout": WIDGETS_LAYOUT
            }
        }
    }
}

# onset_method = "poisson_cusum_bootstrap"
# onset_method_params = {
#     "bg_start": 0,
#     "bg_end": 12,
#     "bootstraps": 1000,
#     "cusum_minutes": 60,
#     "sample_size": 0.75,
#     "limit_averaging": "4 min"
# }

############### Onset Selection ###############
VIEWINGS_HIERARCHY = ["sun", "north", "south", "asun", "omni"]

############### VDA ###############
C = 299_792_458
AU_TO_M_RATIO = 1.495978707e11

M_REST = {
    "protons": 938.27,
    "electrons": 0.511
}