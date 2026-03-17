from datetime import datetime, timezone


class VDA_parameters:

    def __init__(self):
        self.input_type: int = 0
        self.date_start: datetime = datetime(2021, 10, 28, 14, 0)
        self.date_end: datetime = datetime(2021, 10, 28, 20, 0)
        self.date_range_filepath: str = "examples/datetime_range_example.csv"
        self.reference_times_filepath: str = "examples/reference_times_example.csv"
        self.bg_hours_prior: int = 2
        self.bg_hours_after: int = 5
        self.load_data: bool = False
        self.load_data_filepath: str = ""
        self.save_data: bool = False
        self.save_data_filepath: str = ""
        self.viewings_tt: list = [True if v == "sun" else False for v in self.AVAILABLE_VIEWINGS]
        self.resample_frequency: str = "5min"
        self.default_channel_groups: dict = {
            "protons": {
                "HET": [
                    [1, 2, 3],
                    [10, 11, 12],
                    [13, 14, 15],
                    [16, 17, 18],
                    [19, 20, 21],
                    [22, 23, 24],
                    [25, 26, 27],
                    [28, 29, 30, 31]
                ]
            },
            "electrons": {
                "HET": [
                    [0, 1],
                    [2, 3]
                ]
            }
        }
        self.channel_groups: dict = {}
        self.onset_method: str = list(self.AVAILABLE_ONSET_METHODS.keys())[0]
        self.onset_method_parameters: dict = {
            k: v["default"]
            for k, v in self.AVAILABLE_ONSET_METHODS[self.onset_method].items()
        }
        self.onset_selection: int = 0
        self.selected_onsets: dict | None = None
        self.view_dfs: bool = True

    @property
    def sensors(self):
        return set([spec["sensor"] for g in self.channel_groups.values() for spec in g.values()])

    @property
    def particles(self):
        return list(self.channel_groups.keys())

    @property
    def sensors_particles(self):
        sp = {}
        for p, g in self.channel_groups.items():
            for spec in g.values():
                sensor = spec["sensor"]
                try:
                    sp[sensor].append(p)
                except KeyError:
                    sp[sensor] = []
                    sp[sensor].append(p)
                sp[sensor] = list(set(sp[sensor]))
        return sp

    @property
    def viewings(self):
        return [v for i, v in enumerate(self.AVAILABLE_VIEWINGS) if self.viewings_tt[i]]

    @property
    def AVAILABLE_SENSORS(self):
        return ["het", "ept"]

    @property
    def AVAILABLE_PARTICLES(self):
        return ["protons", "electrons"]

    @property
    def AVAILABLE_SENSORS_PARTICLES(self):
        return {
            "het": ["protons", "electrons"],
            "ept": ["protons", "electrons"]
        }

    @property
    def AVAILABLE_CHANNELS(self):
        return {
            "het": {
                "protons": list(range(36)),
                "electrons": list(range(4))
            },
            "ept": {
                "protons": list(range(64)),
                "electrons": list(range(34))
            }
        }

    @property
    def AVAILABLE_VIEWINGS(self):
        return ["sun", "asun", "north", "south", "omni"]

    @property
    def AVAILABLE_ONSET_METHODS(self):
        return {
            "sigma": {
                "s": {
                    "type": int,
                    "min": 1,
                    "max": 5,
                    "default": 3,
                    "description": "Threshold (<this parameter>*<standard deviation>):",
                },
                "n": {
                    "type": int,
                    "min": 1,
                    "max": 5,
                    "default": 3,
                    "description": "Number of consecutive points that should cross the threshold:",
                },
                "bg_start": {
                    "type": int,
                    "min": 0,
                    "max": 9999,
                    "default": 0,
                    "description": "Point index to start the background sampling:",
                },
                "bg_end": {
                    "type": int,
                    "min": 1,
                    "max": 10000,
                    "default": 12,
                    "description": "Point index to end the background sampling:",
                },
            }
        }

# ############### Particle Data ###############

# # onset_method = "poisson_cusum_bootstrap"
# # onset_method_params = {
# #     "bg_start": 0,
# #     "bg_end": 12,
# #     "bootstraps": 1000,
# #     "cusum_minutes": 60,
# #     "sample_size": 0.75,
# #     "limit_averaging": "4 min"
# # }
