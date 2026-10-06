from dataclasses import dataclass, field, fields
from datetime import datetime
from enum import IntEnum
from typing import ClassVar

import pandas as pd


class OnsetSelection(IntEnum):
    USE_ALL = 0
    INTERACTIVE = 1


AVAILABLE_SENSORS_PARTICLES = {
    "het": ("protons", "electrons"),
    "ept": ("protons", "electrons"),
}

AVAILABLE_CHANNELS = {
    "het": {
        "protons": tuple(range(36)),
        "electrons": tuple(range(4)),
    },
    "ept": {
        "protons": tuple(range(64)),
        "electrons": tuple(range(34)),
    },
}

AVAILABLE_VIEWINGS = ("sun", "asun", "north", "south", "omni")

AVAILABLE_ONSET_METHODS = {
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
    }
}

# Parameters removed in v0.3.0 and what replaces them
REMOVED_PARAMETERS = {
    "input_type": "Set input_filepath to the events file (its type is deduced from its columns), or leave it empty to use date_ranges.",
    "date_range_filepath": "Set input_filepath to the events file.",
    "reference_times_filepath": "Set input_filepath to the events file.",
    "date_start": "Set date_ranges to the list of (start, end) datetime ranges of the events.",
    "date_end": "Set date_ranges to the list of (start, end) datetime ranges of the events.",
    "viewings_tt": 'Set viewings to a list of viewing names, e.g. ["sun", "north"].',
}


def _default_onset_method_parameters(method: str) -> dict:
    return {k: v["default"] for k, v in AVAILABLE_ONSET_METHODS[method].items()}


@dataclass
class VDA_parameters:
    # File with the events (datetime ranges or reference times). If empty, date_ranges are used
    input_filepath: str = ""
    # (start, end) datetime range of each event, used when there is no events file
    date_ranges: list = field(default_factory=lambda: [(datetime(2021, 10, 28, 14, 0), datetime(2021, 10, 28, 20, 0))])
    # Data range of the files with reference times, in hours before and after the reference time
    bg_hours_prior: int = 2
    bg_hours_after: int = 5
    # Default background window of the events without one in the events file, in minutes after the start time
    bg_after_start: tuple = (0, 60)
    load_data: bool = False
    load_data_filepath: str = ""
    save_data: bool = False
    save_data_filepath: str = ""
    # Selected viewings. Their order is the priority of the "Use all" onset selection
    viewings: list = field(default_factory=lambda: ["sun"])
    resample_frequency: str = "5min"
    default_channel_groups: dict = field(default_factory=lambda: {
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
    })
    channel_groups: dict = field(default_factory=dict)
    onset_method: str = next(iter(AVAILABLE_ONSET_METHODS))
    onset_method_parameters: dict = field(
        default_factory=lambda: _default_onset_method_parameters(next(iter(AVAILABLE_ONSET_METHODS)))
    )
    onset_selection: OnsetSelection = OnsetSelection.USE_ALL
    # Selected viewing of each grouped channel, set by the onset selection
    selected_onsets: pd.DataFrame | None = field(default=None, repr=False, compare=False)
    view_dfs: bool = True

    AVAILABLE_SENSORS_PARTICLES: ClassVar[dict] = AVAILABLE_SENSORS_PARTICLES
    AVAILABLE_CHANNELS: ClassVar[dict] = AVAILABLE_CHANNELS
    AVAILABLE_VIEWINGS: ClassVar[tuple] = AVAILABLE_VIEWINGS
    AVAILABLE_ONSET_METHODS: ClassVar[dict] = AVAILABLE_ONSET_METHODS

    def __setattr__(self, name, value):
        # catches removed parameters and typos, which would otherwise be silently ignored
        if name in REMOVED_PARAMETERS:
            raise AttributeError(f"The {name} parameter was removed in v0.3.0. {REMOVED_PARAMETERS[name]}")
        if name not in {f.name for f in fields(self)}:
            raise AttributeError(f"VDA_parameters has no parameter '{name}'")
        if name == "onset_selection":
            try:
                value = OnsetSelection(value)
            except ValueError:
                options = ", ".join(f"OnsetSelection.{s.name} ({s.value})" for s in OnsetSelection)
                raise ValueError(f"Unknown onset selection {value!r}. Use one of: {options}") from None
        super().__setattr__(name, value)

    @property
    def sensors_particles(self):
        sp = {}
        for p, g in self.channel_groups.items():
            for spec in g.values():
                particles = sp.setdefault(spec["sensor"], [])
                if p not in particles:
                    particles.append(p)
        return sp
