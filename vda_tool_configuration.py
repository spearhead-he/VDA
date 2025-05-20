from datetime import datetime, timezone


class VDA_parameters:

    def __init__(self):
        self.input_type: int = 0
        self.date_start: datetime = datetime(2021, 5, 22, 19, 45)
        self.date_end: datetime = datetime(2021, 5, 23, 2, 45)
        self.date_range_filepath: str = "examples/datetime_range_example.csv"
        self.reference_times_filepath: str = "examples/reference_times_example.csv"
        self.bg_hours_prior: int = 2
        self.bg_hours_after: int = 5
        self.load_data: bool = False
        self.load_data_filepath: str = ""
        self.save_data: bool = False
        self.save_data_filepath: str = ""
        self.sensors_tt: list = [True for _ in self.AVAILABLE_SENSORS]
        self.particles_tt: list = [True for _ in self.AVAILABLE_PARTICLES]
        self.sensors_particles_tt: dict = {
            s: [True for _ in self.AVAILABLE_SENSORS_PARTICLES[s]]
            for s in self.AVAILABLE_SENSORS_PARTICLES.keys()
        }
        self.viewings_tt: list = [True for _ in self.AVAILABLE_VIEWINGS]
        self.resample_frequency: str = ""
        self.group_sizes: dict = {
            s: {
                p: 2
                for p in ps
            } 
            for s, ps in self.AVAILABLE_SENSORS_PARTICLES.items()
        }
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
        return [s for i, s in enumerate(self.AVAILABLE_SENSORS) if self.sensors_tt[i]]

    @property
    def particles(self):
        return [
            p for i, p in enumerate(self.AVAILABLE_PARTICLES) if self.particles_tt[i]
        ]

    @property
    def sensors_particles(self):
        return {
            s: [p for i, p in enumerate(self.AVAILABLE_SENSORS_PARTICLES[s]) if tt[i]]
            for s, tt in self.sensors_particles_tt.items()
        }

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
