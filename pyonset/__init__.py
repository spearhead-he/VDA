# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    pass  # package is not installed

import datetime
import glob
import os
import pickle
import warnings

import astropy.units as u
import ipywidgets as widgets
import matplotlib.colors as cl
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from pandas.tseries.frequencies import to_offset
from astropy import constants as const
# from matplotlib.ticker import LogLocator, ScalarFormatter
from matplotlib import dates as mdates
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.dates import DateFormatter
from matplotlib.offsetbox import AnchoredText
from sunpy.util.net import download_file

from sunpy import __version__

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

from IPython.display import Markdown

import seppy.tools # this is imported first to avert ImportError in the data loaders that are imported after
import seppy.util as util
from seppy.loader.psp import (calc_av_en_flux_PSP_EPIHI,
                              calc_av_en_flux_PSP_EPILO, psp_isois_load)
from seppy.loader.soho import calc_av_en_flux_ERNE, soho_load
from seppy.loader.solo import epd_load
from seppy.loader.stereo import calc_av_en_flux_HET as calc_av_en_flux_ST_HET
from seppy.loader.stereo import calc_av_en_flux_SEPT, stereo_load
from seppy.loader.wind import wind3dp_load
from seppy.tools import Event

"""
A library that holds the Onset, BackgroundWindow and OnsetStatsArray classes.

@Author: Christian Palmroos <chospa@utu.fi>

@Updated: 2024-12-04

Known problems/bugs:
    > Does not work with SolO/STEP due to electron and proton channels not defined in all_channels() -method
"""


__author__ = "Christian Palmroos"
__email__ = "chospa@utu.fi"

# Some useful global constants
CURRENT_PATH = os.getcwd()
STANDARD_FIGSIZE = (21,9)
VDA_FIGSIZE = (16,9)
C_SQUARED = const.c.value*const.c.value
BACKGROUND_ALPHA = 0.15 # used for the background shading when plotting
TITLE_FONTSIZE = 30
AXLABEL_FONTSIZE = 26
TICK_LABELSIZE = 22
TXTBOX_SIZE = 23
LEGEND_SIZE = 24

COLOR_SCHEME = {
    "median" : "red",
    "mode" : "navy",
    "mean" : "darkorange",
    "1-sigma" : "red",
    "2-sigma" : "blue"
}

ELECTRON_IDENTIFIERS = ("electrons", "electron", 'e')
PROTON_IDENTIFIERS = ("protons", "proton", "ions", "ion", 'p', 'i', 'H')

SEPPY_SPACECRAFT = ("sta", "stb", "solo", "psp", "wind", "soho")
SEPPY_SENSORS = {"sta" : ("sept", "het"),
                 "stb" : ("sept", "het"),
                 "solo" : ("ept", "het"),
                 "psp" : ("isois_epilo", "isois_epihi"),
                 "wind" : ("3dp"),
                 "soho" : ("erne-hed", "ephin")
                 }

# We recommend to have at least this many data points in the background for good statistics
MIN_RECOMMENDED_POINTS = 100

# SOHO / EPHIN e300 channel is invalid from this date onwards
EPHIN_300_INVALID_ONWARDS = pd.to_datetime("2017-10-04 00:00:00")

CUSUM_WINDOW_RESOLUTION_MULTIPLIERS = (4,8,16,32)

# A new class to inherit everything from serpentine Event, and to extend its scope of functionality
class Onset(Event):

    def __init__(self, start_date, end_date, spacecraft, sensor, species, data_level, data_path, viewing=None, radio_spacecraft=None, threshold=None,
                 data=None, unit=None):

        # By default we download data, not provide it
        if data is None:
            super().__init__(start_date, end_date, spacecraft, sensor,
                    species, data_level, data_path, viewing, radio_spacecraft,
                    threshold)
            self.custom_data = False
            self.unit = r"Intensity [1/(cm$^{2}$ sr s MeV)]" if unit is None else unit

            # Check here that the spacecraft and instrument are SEPpy-compatible.
            # to provide the custom data.
            if self.spacecraft not in SEPPY_SPACECRAFT:
                raise Exception(f"Note that {self.spacecraft} is not a SEPpy-compatible spacecraft. You need to provide your own data with the keyword: 'data'.")

            # if self.sensor not in SEPPY_SENSORS[self.spacecraft]:
            #     raise Exception(f"Note that {self.sensor} is not a recognized instrument of {self.spacecraft}. You need to provide your own data with the keyword: 'data'.")

        else:

            self.start_date = start_date
            self.end_date = end_date
            self.spacecraft = spacecraft
            self.sensor = sensor
            self.species = species
            self.data_level = data_level
            self.data_path = data_path
            self.viewing = viewing
            self.radio_spacecraft = radio_spacecraft
            self.threshold = threshold

            if not isinstance(data,pd.DataFrame):
                raise TypeError("Custom data needs to be provided as a DataFrame indexed by time!")
            self.data = data.copy(deep=True)
            self.current_df_e = self.data
            self.current_df_i = self.current_df_e
            self.unit = r"Intensity [1/(cm$^{2}$ sr s MeV)]" if unit is None else unit

            # Custom data flag prevents SEPpy functions from being called, as they would cause errors
            self.custom_data = True

            # Lets the user know that the object is initialized with custom settings
            print("Utilizing user-input data. Some SEPpy functionality may not work as intended.")

        # Everytime an onset is found any way, the last used channel should be updated
        self.last_used_channel = np.nan

        # The background window is stored to this attribute when cusum_onset() is called with a BootStrapWindow input
        self.background = np.nan

        # This list is for holding multiple background windows if such were to be used
        self.list_of_bootstrap_windows = []
        self.window_colors = ["blue", "orange", "green", "purple", "navy", "maroon"]

        # This will turn true once the extensive statistics analysis is run
        self.mean_of_medians_onset_acquired = False

        # Choosing the particle species identifier for the titles etc
        if self.species in ["electron", "electrons", 'e']:
            self.species = 'e'
            self.s_identifier = "electrons"
        elif self.species in ["proton", "protons", 'p', 'H']:
            self.species = 'p'
            self.s_identifier = "protons"
        elif self.species in ["ion", "ions", 'i']:
            self.species = 'p'
            self.s_identifier = "ions"
        else:
            self.s_identifier = self.species

        # This is a check to make sure viewing is correctly set for single-aperture instruments
        if self.viewing:
            self.check_viewing()

        self.all_channels = {
            "sta_let_e" : np.arange(1, dtype=int),
            "sta_sept_e" : np.arange(2,17, dtype=int),
            "sta_het_e" : np.arange(3, dtype=int),
            "sta_let_p" : np.arange(1, dtype=int),
            "sta_sept_p" : np.arange(2,32, dtype=int),
            "sta_het_p" : np.arange(11, dtype=int),
            
            "stb_let_e" : np.arange(1, dtype=int),
            "stb_sept_e" : np.arange(2,17, dtype=int),
            "stb_het_e" : np.arange(3, dtype=int),
            "stb_let_p" : np.arange(1, dtype=int),
            "stb_sept_p" : np.arange(2,32, dtype=int),
            "stb_het_p" : np.arange(11, dtype=int),

            "solo_step_e" : np.arange(1, dtype=int),
            "solo_ept_e" : np.arange(34, dtype=int),
            "solo_het_e" : np.arange(4, dtype=int),
            "solo_step_p" : np.arange(1, dtype=int),
            "solo_ept_p" : np.arange(64, dtype=int),
            "solo_het_p" : np.arange(36, dtype=int),

            "soho_ephin_e" : (150, 300, 1300, 3000),
            "soho_erne_p" : np.arange(10, dtype=int),

            "wind_3dp_e" : np.arange(7, dtype=int),
            "wind_3dp_p" : np.arange(9, dtype=int),

            "psp_isois-epihi_e" : np.arange(19, dtype=int),
            "psp_isois-epihi_p" : np.arange(15, dtype=int),

            "psp_isois-epilo_e" : np.arange(12, dtype=int),
        }

        # This dictionary holds the median onset time and weighted confidence intervals derived by CUSUM-bootstrap 
        # hybrid method. In the dictionary each channel identifier is paired with a list of [mode, median, 1st_minus, 1st_plus, 2nd_minus, 2nd_plus]
        self.onset_statistics = {}

        # This is a dictionary that holds the information of each instrument's minimum time resolution. That value
        # will be used in the case that cusum-bootstrap yields uncertainty smaller than that.
        self.minimum_cadences = {
            "bepicolombo_sixs-p" : pd.Timedelta("8 s"),
            "psp_isois-epihi" : pd.Timedelta("1 min"),
            "psp_isois-epilo" : pd.Timedelta("1 min"),
            "soho_erne" : pd.Timedelta("1 min"),
            "soho_ephin" : pd.Timedelta("1 min"),
            "sta_sept" : pd.Timedelta("1 min"),
            "sta_het" : pd.Timedelta("1 min"),
            "stb_sept" : pd.Timedelta("1 min"),
            "stb_het" : pd.Timedelta("1 min"),
            "solo_step" : pd.Timedelta("1 s"),
            "solo_ept" : pd.Timedelta("1 s"),
            "solo_het" : pd.Timedelta("1 s"),
            "wind_3dp" : pd.Timedelta("12 s")
        }

    def check_viewing(self, returns=False):
        """
        A method that checks at initialization that a viewing direction wasn't erroneously appointed
        for an instrument that doesn't have viewing directions.
        Sets self.viewing to None for STA/STB HET and all SOHO instruments (ERNE and EPHIN)

        returns : {bool}
                    Switch to return the viewing
        """

        if self.sensor.lower() in ("het", "erne", "ephin"):

            # STA/STB HET, and SOHO all instruments don't have viewings -> assign None to them
            if self.spacecraft.lower() in ("sta", "stb", "soho"):
                self.viewing = None
            else:
                self.viewing = self.viewing

        if returns:
            return self.viewing

    def __repr__(self):
        return str(f"({self.spacecraft},{self.sensor},{self.species})")

    def get_all_channels(self):
        """
        Returns a range(first,last+1) of all the channel identifiers for any unique instrument+species pair.

        KeyError is thrown if the spacecraft+sensor+species -combination doesn't correspond to any of the
        recognized channel lists. This is handled by using the self.data (pd.DataFrame) column names.
        """

        try:
            return self.all_channels[f"{self.spacecraft}_{self.sensor}_{self.species}"]
        except KeyError:
            return np.array(self.data.columns)

    def get_current_channel_str(self):
        """
        Gets the string for the current energy channel
        """
        return self.get_channel_energy_values(returns="str")[self.last_used_channel]

    def get_minimum_cadence(self):
        try:
            return self.minimum_cadences[f"{self.spacecraft}_{self.sensor}"]
        except KeyError:
            # return pd.Timedelta(self.data.index.freq)
            return pd.Timedelta(get_time_reso(self.data))

    def get_time_resolution_str(self, resample):
        # Choose resample as the averaging string if it exists
        if resample:
            time_reso_str = f"{resample} averaging" 
        # If not, check if the native cadence is less than 60. If yes, address the cadence in seconds
        elif self.get_minimum_cadence().seconds<60:
            time_reso_str = f"{self.get_minimum_cadence().seconds} s data"
        # If at least 60 seconds, address the cadence in minutes
        else:
            try:
                time_reso_str = f"{int(self.get_minimum_cadence().seconds/60)} min data"
            except ValueError:
                time_reso_str = "Unidentified time resolution"

        return time_reso_str

    def get_custom_channel_energies(self):
        """
        Gets the energy values that have been set by the user.
        Returns energy channel low and high boundaries in eVs.

        Handles the AttributeError that is caused by non-defined energy boundary
        values by printing the error, an explanation of how to fix it and returning None.
        """

        try:
            return self.channel_energy_lows, self.channel_energy_highs
        except AttributeError as e:
            print(e)
            print("This is caused by user-defined custom data, that has no defined energy channel boundary values.")
            print("Define energy channel boundaries with 'set_custom_channel_energies()'!")
            return None

    def add_bootstrap_window(self, window):
        """
        A little helper method to add bootstrapping windows to the class list
        """
        self.list_of_bootstrap_windows.append(window)

    def update_onset_statistics(self, stats):
        """
        Adds the weighted median and confidence intervals of a specific channel to onset_statistics dictionary
        """

        self.onset_statistics[self.last_used_channel] = stats

    def input_nat_onset_stats(self, key):
        """
        This method can be used to e.g., fill in SolO/EPT channel 0 stats, which is unavailable in latest data.
        It may also be used to erase clearly erroneous onset time data from the object.
        """
        nan_stats = [pd.NaT, pd.NaT, pd.NaT, pd.NaT, pd.NaT, pd.NaT]
        self.onset_statistics[key] = nan_stats

    def set_data_frequency(self, freq):
        """
        Tries to set a frequency for the data. Useful especially incase of custom data.
        
        Parameters:
        -----------
        freq : {str}
                A pandas-compatible time string that represents the frequency of 
                the time series data, e.g., '24 s' or '1 min'.
        """
        self.data.index.freq = freq

    def set_channel_strings(self, channel_lows, channel_highs, unit):
        """
        Sets the dictionary of channel energy strings, e.g., '45 - 55 keV'.
        """
        channel_en_dict = {}
        for i, low in enumerate(channel_lows):
            
            channel_en_dict[self.data.columns[i]] = f"{low} - {channel_highs[i]} {unit}"
        
        self.channel_en_dict = channel_en_dict

    def set_custom_channel_energies(self, low_bounds, high_bounds, unit="MeV"):
        """
        Sets the channel energy boundary values. Automatically updates self.channel_en_dict
        
        Parameters:
        -----------
        low_bounds : {array-like}
        
        high_bounds : {array_like}
        
        unit : {str} default 'MeV'. Choose either 'eV', 'keV' or 'MeV'.
        """

        if len(low_bounds) != len(self.data.columns) or len(high_bounds) != len(self.data.columns):
            raise ValueError("Either low_bounds or high_bounds has an incorrect amount of entries!")

        if unit not in ("eV", "keV", "MeV"):
            raise ValueError(f"Unit {unit} does not appear to be any of the recognized units of energy! Choose either 'eV', 'keV' or 'MeV'.")

        if unit == "MeV":
            self.channel_energy_lows = np.array([low*u.MeV.to("eV") for low in low_bounds])
            self.channel_energy_highs = np.array([high*u.MeV.to("eV") for high in high_bounds])

        if unit == "keV":
            self.channel_energy_lows = np.array([low*u.keV.to("eV") for low in low_bounds])
            self.channel_energy_highs = np.array([high*u.keV.to("eV") for high in high_bounds])

        if unit == "eV":
            self.channel_energy_lows = np.array([low for low in low_bounds])
            self.channel_energy_highs = np.array([high for high in high_bounds])
        
        self.set_channel_strings(channel_lows=low_bounds, channel_highs=high_bounds, unit=unit)


    def cusum_onset(self, channels, background_range, viewing=None, resample=None, cusum_minutes=30, sigma=2, title=None, save=False, savepath=None, 
                    yscale='log', ylim=None, erase=None, xlim=None, show_stats=True, diagnostics=False, plot=True, fname:str=None):
        """
        Does a Poisson-CUSUM-method-based onset analysis for given OnsetAnalysis object.
        Based on an earlier version by: Eleanna Asvestari <eleanna.asvestari@helsinki.fi>

        Parameters
        ----------
        channel: int or a list of 2 ints
                Index of the channel inside the DataFrame that holds the data.
        background_range : tuple/list or {BootstrapWindow}
                The starting and ending points of what is considered background. 
        resample: str
                Pandas-compatible time string, e.g., '1min' or '30s'
        cusum_minutes: int
                The amount of minutes the intensity should stay above threshold until onset is identified. Corresponds to the amount
                of consecutive "out-of-control" alarms one should get before finding the onset.
        sigma: int, default 2
                How many standard deviations is the mu_d variable of the CUSUM function?
        title: str, default None
                A custom title for your plot.
        save: bool
                Wether or not to save the figure in local "plots" directory. Only works if plot=True.
        savepath: str, default None
                The path to which to save the plot, if not the default "plots".
        yscale: str
                The scale of the y-axis, either 'log' or 'linear'.
        ylim: 2-tuple, default None
                Defines the lower and upper bounds for the plots vertical axis.
        erase: 2-tuple, default None
                If there are spikes in the background one wishes to omit, set the threshold value and ending point for ignoring those
                values here. Example: [1000, '2022-05-10 12:00'] discards all values 1000 or above until 2022-05-10 noon.
        xlim: 2-tuple (str, str), default None
                Pandas-compatible dtateime strings. Defines the start and end for the horizontal axis of the plot.
        show_stats : bool, default True
                Switch to show the mean of medians onset and confidence intervals. Also shows the earliest possible onset (first nonzero CUSUM value)
        diagnostics: bool, default False
                Diagnostic tool allows for plotting the CUSUM function, and printing many 'behind the curtains' variables for the user to see.
        plot: bool, default True
                Switch to produce a plot. 
        fname : {str} default None
                A custom name for the figure if saved.

        Returns:
        ---------
        onset_stats: list
                A list of onset statistics. [background_mean, background_mean+sigma*std, k_round, normed_values, CUSUM_function, Timestamp_of_Onset]
        flux_series: Pandas Series
                The series used for producing the plot and running the onset analysis.
        """

        spacecraft = self.spacecraft.upper()

        # Set the viewing even in the case it was not input by a user
        viewing = viewing if viewing else self.viewing

        # Logic of this check:
        # If viewing was set (meaning viewing evaluates to True) and Onset.check_viewing() returns None 
        # (which evaluates to False), then viewing must be invalid.
        if viewing and not self.check_viewing(returns=True):
            raise ValueError("Invalid viewing direction!")

        color_dict = {
            'onset_time': '#e41a1c',
            'bg_mean': '#e41a1c',
            'flux_peak': '#1a1682',
            'bg': '#de8585'
        }

        # This check was initially transforming the 'channels' integer to a tuple of len==1, but that
        # raised a ValueError with solo/ept. However, a list of len==1 is somehow okay.
        # Also save the channel to a class attribute for later use
        if isinstance(channels, (int, np.int64)):
            self.last_used_channel = int(channels)
            channels = [int(channels)]

        if isinstance(background_range,(list,tuple)):
            background_start, background_end = pd.to_datetime(background_range[0]), pd.to_datetime(background_range[1])
        
        if isinstance(background_range, BootstrapWindow):
            background_start, background_end = background_range.start, background_range.end
            self.background = background_range

        # By default we do not use custom data
        if not self.custom_data:
            flux_series, en_channel_string = self.choose_flux_series(channels=channels, viewing=viewing)
        else:
            flux_series, en_channel_string = self.data[channels], channels
            self.last_used_channel = channels

        # Glitches from the data should really be erased BEFORE resampling data
        if erase is not None:
            flux_series, glitches = erase_glitches(flux_series, erase[0], erase[1])
        else:
            glitches = None

        # Resample data if requested
        if resample is not None:
            flux_series = util.resample_df(flux_series, resample)

        # Here just to make sure, check that the flux series is cut correctly
        if len(flux_series) == 0:
            raise Exception("Intensity data selection is empty! Check that plot limits are correct.")

        # The series is indexed by time
        time = flux_series.index

        # Saves the intensity time series to a class attribute, handy if debugging or such is needed
        self.flux_series = flux_series

        # Onset is not yet found
        onset_found = False

        # Background stats are a tuple of (mean, standard_deviation)
        background_stats = calculate_mean_and_std(background_start, background_end, flux_series)

        # Also check that the background parameters are calculated using a reasonable background_window
        if background_stats[0] == 0 or np.isnan(background_stats[0]):
            background_warning = f"WARNING: background mean is {background_stats[0]}! If this was unintended, it may have been caused by an incorrect choice of background window."
            warnings.warn(message=background_warning)

        # Calculate the length of cusum_window in terms of data points from the time resolution
        time_reso = get_time_reso(flux_series)
        cusum_window = calculate_cusum_window(time_reso, cusum_minutes)

        # The function finds the onset and returns a list of stats related to the onset
        # onset_stats = [ma, md, k_round, h, norm_channel, cusum, onset_time]
        if self.unit not in ("Counting rate [1/s]", "Count rate [1/s]", "Counting rate", "Count rate", "counting rate", \
                             "count rate", "Counting_rate", "Count_rate", "Counts", "counts"):
            onset_stats = onset_determination(background_stats, flux_series, cusum_window, background_end, sigma)
        else:
            # If the unit is count rate (1/s), then employ Poisson-CUSUM without using z-standardized intensity
            onset_stats = onset_determination_cr(background_stats, flux_series, cusum_window, background_end, sigma)

        # If the timestamp of onset is not NaT, then onset was found
        if not isinstance(onset_stats[-1],pd._libs.tslibs.nattype.NaTType):
            onset_found = True

        # Prints out useful information if diagnostics is enabled
        if diagnostics:
            print(f"Cusum_window, {cusum_minutes} minutes = {cusum_window} datapoints")
            print(f"onset time: {onset_stats[-1]}")
            print(f"mu and sigma of background intensity: \n{np.round(background_stats[0],2)}, {np.round(background_stats[1],2)}")

        # --Only plotting related code from this line onward ->

        # Before starting the plot, save the original rcParam options and update to new ones
        original_rcparams = self.save_and_update_rcparams("onset_tool")

        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        # Setting the x-axis limits
        if xlim is None:
            xlim = [time[0], time[-1]]
        else:

            # Check that xlim makes sense
            if xlim[0] == xlim[1]:
                raise Exception("xlim[0] and xlim[1] cannot be the same time!")

            xlim[0], xlim[1] = pd.to_datetime(xlim[0]), pd.to_datetime(xlim[1])

        # Just make sure that the user did not mess up their x-axis limits and the background for the cusum function
        if (background_range is not None) and (xlim is not None):
            # Check if background is separated from plot range by over a day, issue a warning if so, but don't
            if (background_start < xlim[0] - datetime.timedelta(days=1) and background_start < xlim[1] - datetime.timedelta(days=1)) or \
               (background_end > xlim[0] + datetime.timedelta(days=1) and background_end > xlim[1] + datetime.timedelta(days=1)):
                background_warning = "NOTICE that your background_range is separated from plot_range by over a day.\nIf this was intentional you may ignore this warning."
                warnings.warn(message=background_warning)

        # Figure limits and scale
        flux_in_plot = flux_series.loc[(flux_series.index < xlim[-1])&(flux_series.index >= xlim[0])]
        ylim = set_fig_ylimits(ax=ax, flux_series=flux_in_plot, ylim=ylim)
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_yscale(yscale)

        # The measurements. kw where dictates where the step happens -> "mid" for at the middle of the bin
        ax.step(time, flux_series.values, c="C0", where="mid")

        if erase is not None:
            ax.scatter(glitches.index, glitches.values, s=3, c='maroon')
            if diagnostics:
                print("omitted values:")
                print(glitches)

        # These are for bughunting
        if diagnostics:
            ax.step(time, onset_stats[-3], color="darkgreen", label=r"$I_{z-score}$")
            ax.step(time, onset_stats[-2], color="maroon", label="CUSUM")
            ax.axhline(y=onset_stats[2], ls="--", color='k', label='k')
            ax.axhline(y=onset_stats[3], ls="-.", color='k', label='h')

        # Onset time
        if onset_found:
            ax.axvline(x=onset_stats[-1], linewidth=1.5, color=color_dict["onset_time"], linestyle= '-', 
                        label="onset time")

            # Textbox indicating onset time
            onset_timestr = onset_stats[-1].strftime("%H:%M:%S")
            plabel = AnchoredText(f"Onset time: {onset_timestr} ", prop=dict(size=TXTBOX_SIZE), frameon=True, loc="lower right") #str(onset_stats[-1])[:19]
            plabel.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            plabel.patch.set_linewidth(2.0)
            ax.add_artist(plabel)

        # Background mean
        ax.axhline(y=onset_stats[0], linewidth=2, color=color_dict["bg_mean"], linestyle= "--", 
                    label=r'$\mu_{a}$ for background intensity')

        # Background mean + n*std
        ax.axhline(y=onset_stats[1], linewidth=2, color=color_dict["bg_mean"], linestyle= ':',
                    label=r'$\mu_{d}$ for background intensity')

        # Background shaded area
        if isinstance(background_range, BootstrapWindow):
            background_range.draw_background(ax=ax)
        else:
            ax.axvspan(xmin=background_start, xmax=background_end,
                        color=color_dict["bg"], label="Background", alpha=0.15)

        # Plot the earliest possible onset time, mean of medians and confidence intervals if they're found
        if show_stats:

            # Mean of median onsets accompanied by their confidence intervals
            if self.mean_of_medians_onset_acquired:

                ax.axvline(x=self.mean_of_medians_onset, linewidth=1.5, color="magenta", linestyle='-',
                        label="mean of medians onset")

                ax.axvspan(xmin=self.conf1_low, xmax=self.conf1_high, color="red", alpha=0.3, label="~68 % confidence")
                ax.axvspan(xmin=self.conf2_low, xmax=self.conf2_high, color="blue", alpha=0.3, label="~95 % confidence")


        ax.set_xlabel(f"Time ({time[0].year})", fontsize=AXLABEL_FONTSIZE)
        ax.set_ylabel(self.unit, fontsize=AXLABEL_FONTSIZE)

        # Date tick locator and formatter
        ax.xaxis_date()
        set_standard_ticks(ax=ax)

        utc_dt_format1 = DateFormatter('%H:%M \n%m-%d')
        ax.xaxis.set_major_formatter(utc_dt_format1)

        # Setting the title
        if title is None:

            if viewing:
                ax.set_title(f"{spacecraft}/{self.sensor.upper()} ({viewing}) {en_channel_string} {self.s_identifier}\n{time_reso} data", fontsize=TITLE_FONTSIZE)
            else:
                ax.set_title(f"{spacecraft}/{self.sensor.upper()} {en_channel_string} {self.s_identifier}\n{time_reso} data", fontsize=TITLE_FONTSIZE)

        else:
            ax.set_title(title, fontsize=TITLE_FONTSIZE)

        if diagnostics:
            ax.legend(loc="best", fontsize=LEGEND_SIZE)

        # Attach the figure to class attribute even if not saving the figure
        self.fig, self.ax = fig, ax

        if save:
            if savepath is None:
                savepath = CURRENT_PATH

            # A custom name for the figure
            if fname is not None:
                fig.savefig(fname=f"{savepath}{os.sep}{fname}",
                            facecolor="white", transparent=False, bbox_inches="tight")

            # Use a default name for the figure
            else:

                if spacecraft.lower() in ["bepicolombo", "bepi"]:
                    plt.savefig(f"{savepath}{os.sep}{self.spacecraft}{self.sensor}_side{viewing}_{self.species}_{channels}_onset.png", transparent=False,
                            facecolor='white', bbox_inches='tight')
                elif viewing != "" and viewing is not None:
                    plt.savefig(f"{savepath}{os.sep}{self.spacecraft}{self.sensor}_{viewing.lower()}_{self.species}_{channels}_onset.png", transparent=False,
                            facecolor='white', bbox_inches='tight')
                else:
                    plt.savefig(f"{savepath}{os.sep}{self.spacecraft}{self.sensor}_{self.species}_{channels[:]}_onset.png", transparent=False,
                            facecolor='white', bbox_inches='tight')

        if plot:
            plt.show()
        else:
            plt.close()

        return onset_stats, flux_series


    def plot_all_channels(self, viewing:str=None, resample:str=None, omit:list=None, cmap:str="twilight_shifted", xlim=None, title:str=None,
                          save=False, savepath=None):
        """
        Plots by default all the energy channels of the currently chosen instrument in the same time series plot.
        
        Parameters:
        -----------
        viewing : {str}, optional
                    The viewing direction of the instrument. If not provided then use the currently 
                    active viewing direction from self.viewing.
        resample : {str}, optional
                    The time-averaging for the time series. Use a pandas-compatible time string, e.g., '1min' or '600s'
        omit : {list, tuple}, optional
                    A list or a tuple of the channels to forgo plotting.
        cmap : {str}, optional
                    The colormap used to plot the channels.
        xlim : {tuple or list with __len__ == 2}, optional
                    The x-axis limits of the figure.
        title : {str}, optional
                    A custom title for the figure.
        save : {bool}, optional
                    Boolean save switch, default False.
        savepath : {str}, optional
                    Either relative or absolute path leading to the directory where the plot should be saved to.
        
        Returns:
        -----------
        fig, ax : The figure and the axis objects.
        """

        # Choose viewing
        viewing = viewing if viewing else self.viewing

        # Collect all channels to a list
        all_channels = self.get_all_channels()

        # Choose which channels to plot, by default plot everything
        plotted_channels = all_channels if not omit else list(set(all_channels) - set(omit))

        # Get a colormap for the different channels.
        colors = plt.get_cmap(cmap, len(plotted_channels))

        # Initializing the figure and axes
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        # Before starting the plot, save the original rcParam options and update to new ones
        # original_rcparams = self.save_and_update_rcparams("onset_tool")
        # This seemed to sometimes get figure visual parameters messed up, such as fontsize of the title and figure size

        # Plotting loop
        for i, channel in enumerate(plotted_channels):

            # Getting a time series from the object's dataframe
            intensity, en_channel_string = self.choose_flux_series(channels=[int(channel)], viewing=viewing)

            if isinstance(en_channel_string,(list, np.ndarray)):
                en_channel_string = en_channel_string[0]

            # If time-averaging, do it here before plotting
            if resample:
                intensity = util.resample_df(intensity, resample)

            # Plotting a channel. kw where=="mid" to ensure that step happens in the middle of the time bin
            ax.step(intensity.index, intensity.values, where="mid", label=en_channel_string, color=colors(i))


        # Setting the x-axis limits
        if not xlim:
            xlim = [intensity.index[0], intensity.index[-1]]
        else:

            # Check that xlim makes sense
            if xlim[0] == xlim[1]:
                raise Exception("xlim[0] and xlim[1] cannot be the same time!")

            xlim[0], xlim[1] = pd.to_datetime(xlim[0]), pd.to_datetime(xlim[1])
        
        ax.set_xlim(xlim)

        ax.set_yscale("log")

        ax.set_xlabel("Time", fontsize=20)
        ax.set_ylabel(f"Intensity [1/(cm^2 sr s MeV)]", fontsize=20)

        # Tickmarks, their size etc...
        ax.tick_params(which='major', length=6, width=2, labelsize=17)
        ax.tick_params(which='minor', length=5, width=1.4)

        # Date tick locator and formatter
        ax.xaxis_date()
        #ax.xaxis.set_major_locator(ticker.AutoLocator()) # ticker.MaxNLocator(9)
        utc_dt_format1 = DateFormatter('%H:%M \n%Y-%m-%d')
        ax.xaxis.set_major_formatter(utc_dt_format1)

        ax.legend(loc=3, bbox_to_anchor=(1.0, 0.01), prop={'size': 16})

        # Setting the title
        if title is None:

            # Choosing the particle species identifier for the title
            if self.species in ["electron", 'e']:
                s_identifier = 'electrons'
            if self.species in ["proton", 'p', 'H']:
                s_identifier = 'protons'
            if self.species in ["ion", 'i']:
                s_identifier = 'ions'

            if viewing:
                ax.set_title(f"{self.spacecraft.upper()}/{self.sensor.upper()} {s_identifier}\n"
                        f"{resample} averaging, viewing: "
                        f"{viewing.upper()}", fontsize=TITLE_FONTSIZE)
            else:
                ax.set_title(f"{self.spacecraft.upper()}/{self.sensor.upper()} {s_identifier}\n"
                        f"{resample} averaging", fontsize=TITLE_FONTSIZE)
        
        else:
            ax.set_title(title, fontsize=TITLE_FONTSIZE)

        if save:
            if savepath is None:
                savepath = CURRENT_PATH

            if self.spacecraft.lower() in ['bepicolombo', 'bepi']:
                fig.savefig(f"{savepath}{os.sep}{self.spacecraft}_{self.sensor}_side{viewing}_{self.species}_all_channels.png", transparent=False,
                        facecolor='white', bbox_inches='tight')
            elif viewing != "" and viewing:
                fig.savefig(f"{savepath}{os.sep}{self.spacecraft}_{self.sensor}_{viewing.lower()}_{self.species}_all_channels.png", transparent=False,
                        facecolor='white', bbox_inches='tight')
            else:
                fig.savefig(f"{savepath}{os.sep}{self.spacecraft}_{self.sensor}_{self.species}_all_channels.png", transparent=False,
                        facecolor='white', bbox_inches='tight')

        return fig, ax
        


    def statistic_onset(self, channels, Window, viewing:str, resample:str=None, erase:tuple=None, sample_size:float=None, cusum_minutes:int=None, 
                        small_windows:str=None, offset_origins:bool=True, detrend=True, sigma_multiplier=2):
        """
        This method looks at a particular averaging window with length <windowlen>, and draws from it
        points <n_bstraps> times. From these <n_bstraps> different distributions of measurements, it 
        calculates the mean and standard deviation for each distribution. Using the acquired distribution of
        different means and standard deviations, it produces a distribution of expected onset times. Finally
        the method moves the averaging window forward if windows has n_shifts>1.
        
        The procedure described above is repeated <n_windows> times, yielding <n_windows>*<n_bstraps> onset times.
        
        Parameters:
        ------------
        channels : int or list
                    Channel or channels that the onset is searched from
        Window: BootstrapWindow object
                    A BootstrapWindow object that holds information about averaging window parameters
        viewing : str
                    The viewing direction of the instrument
        resample : str
                    A pandas-compatible time string representing the averaging of data
        erase : tuple(float, str)
                    The maximum acceptable intensity value and the timestamp at which larger values are accepted
        sample_size : {float}, optional
                    The proportion of the ensemble of data points inside a window that are drawn for the bootstrapping method
        cusum_minutes : {int}, optional 
                    An integer stating the amount of minutes that threshold exceeding intensity must consistently be observed 
                    before identifying an onset. If not provided, use a set of 4,8, 16 and 32 times the current data resolution
        small_windows : str, default None
                    Choose either 'random' or 'equal' placement of the subwindows. Randomly placed windows can overlap,
                    while equally placed windows start where the previous ended.
        offset_origins : {bool}, default True
                    Boolean switch to control if onset times can be found on a continuous scale or only at intervals allowed by the
                    resolution of the data. By default origin is being offset, i.e., onset can be found at 1-minute resolution
                    regardless of the integration time.
        detrend : bool, default True
                    Controls if the onset times are 'de-trended', i.e., shifted forward in time by 0.5*time_resolution - 30s. Does not 
                    apply to 1min data.
        sigma_multiplier : {int, float} default 2
                    The multiplier for the $\mu_{d}$ variable in the CUSUM method.

        Returns:
        -----------
        mean_onset: a timestamp indicating the mean of all onset times found
        std_onset: the standard deviation of onset times
        most_common_val: the most common onset time found
        percent_mode: the prevalence of the most common onset time (fraction of all onset times)
        onset_list: a list containing all onset times in the order they were calculated by this function
        """
        # Bootstrapping:
        # https://towardsdatascience.com/calculating-confidence-interval-with-bootstrapping-872c657c058d

        # Channels are considered a list in the data loader, that's why this check is done here
        self.last_used_channel = channels
        if isinstance(channels,(int,str)):
            channels = [channels]


        self.current_channel_id = channels[0]

        # Do not search for onset earlier than this point
        big_window_end = Window.end

        # Before going to pick the correct flux_series, check that viewing is reasonable.
        if viewing is not None and not self.check_viewing(returns=True):
            raise ValueError("Invalid viewing direction!")

        # Choose the right intensity time series according to channel and viewing direction.
        # Also remember which channel was examined most recently.
        if not self.custom_data:
            flux_series, self.recently_examined_channel_str = self.choose_flux_series(channels, viewing)
        else:
            flux_series, self.recently_examined_channel_str = self.data[channels[0]], channels[0]

        # Create a timedelta representation of resample, for convenience
        resample_td = pd.Timedelta(resample)
        # By default there will be a list containing timeseries indices of varying origins of offset
        if offset_origins and resample:

            # Separate the unit of resample and its numerical value
            resample_unit_str = find_biggest_nonzero_unit(timedelta=resample_td)
            resample_value = int(resample[:-3]) if resample_unit_str=="min" else int(resample[:-1])
            if resample[-3:] == "min":

                # This integer is one larger than the offset will ever go
                data_res = int(resample_value)

                # There will be a set amount of varying offsets with linear intervals to the time indices of flux_series.
                # offsets are plain integers
                offsets = np.arange(data_res)

                # Using the number of offsets, compile a list of copies of flux_series
                list_of_series = [flux_series if i==0 else flux_series.copy() for i in offsets]

                # Produces a list of indices of varying offsets. The first one in the list has offset=0, which are 
                # just the original indices. These indices will replace the indices of the series in the list.
                list_of_indices = [flux_series.index + to_offset(f"{offset}{resample_unit_str}") for offset in offsets]

                # These are a random numbers to access any of the indices in the list of indices.
                # Calculated beforehand to save computation time.
                random_choices = np.random.randint(data_res, size=Window.bootstraps)

            else:
                list_of_series = [flux_series]
                list_of_indices = [flux_series.index]
                random_choices = np.zeros(Window.bootstraps, dtype=int)

        # If we don't offset origin of timeseries, then just init list_of_series, list_of_indices and random_choises
        # with default values for a single unaltered flux_series
        else:

            list_of_series = [flux_series]
            list_of_indices = [flux_series.index]
            random_choices = np.zeros(Window.bootstraps, dtype=int)


        # Replace the indices of the copies of the series. Each new series will have its indices shifted
        # forward by one minute (in the current version).
        # After the possible shifting of indices, first erase possible erroneous data, only after 
        # resample it if need be.
        for i, series in enumerate(list_of_series):

            # Shifting of indices forward
            series.index = list_of_indices[i]

            # Erasing erroneous values
            if erase:
                series, _ = erase_glitches(series, erase[0], erase[1])

            # Resampling of data to a differrent resolution
            if resample:
                series = util.resample_df(series, resample)

            # Finally placing the series back to the list at the correct position
            list_of_series[i] = series


        # These apply regardless of which series we choose, as they all have the same resolution and hence cusum_window
        if not self.custom_data:
            if not resample:
                time_reso = f"{int(self.get_minimum_cadence().seconds/60)} min" if self.get_minimum_cadence().seconds > 59 else f"{int(self.get_minimum_cadence().seconds)} s"
            else:
                time_reso = get_time_reso(list_of_series[0])

        # in case of custom data: 
        else:
            time_reso = get_time_reso(list_of_series[0])
            if not cusum_minutes:
                raise Exception("Must provide a value for 'cusum_minutes' with custom data!")

        # If cusum_window was not explicitly stated, use a set of multiples of the time resolution as 
        # a set of cusum_windows
        if cusum_minutes:
            cusum_windows = [calculate_cusum_window(time_reso, cusum_minutes)]
        else:
            if time_reso[-3:]=="min":
                cusum_minutes_list = [c*int(time_reso[:-3]) for c in CUSUM_WINDOW_RESOLUTION_MULTIPLIERS]
            else:
                # For now just go with a fixed number of minutes if resolution is less than a minute
                cusum_minutes_list = [c for c in CUSUM_WINDOW_RESOLUTION_MULTIPLIERS]
                # cusum_minutes_list = [int(c*int(time_reso[:-1])/60) for c in cusum_window_resolution_multipliers]
            cusum_windows = [calculate_cusum_window(time_reso, cusum_minutes) for cusum_minutes in cusum_minutes_list]

        # Initialize necessary lists to collect data about onset time and the corresponding background parameters
        onset_list = []
        mus_list = []
        sigmas_list = []

        # Small windows will initialize a set of randomly placed smaller windows inside the chosen (input) 
        # BootstrapWindow object. The small windows will all have n_shifts=0, and bootstraps equal to the
        # amount the original Window had.
        if small_windows:

            # Before doing anything, check that small_windows is the type it is meant to be
            if small_windows not in ("random", "equal"):
                raise Exception(f"small_windows = {small_windows} is not a valid input! Choose either 'random' or 'equal'.")

            # First check that there are shifts in the Window, otherwise it makes no sense to even
            # call subwindows()
            if Window.n_shifts < 1:
                raise Exception(f"The number of small windows can't be {Window.n_shifts}")

            # Produce a list of little windows to loop over and apply bootstrapping
            windowlist = subwindows(Window=Window, window_placement=small_windows)

            # Set the class-attached list to an empty list, in order to stop new windows accumulating there every time
            # this method is run. This is essentially flushing the list clean.
            self.list_of_bootstrap_windows = []
            self.subwindows_plotted = False

        # If no small windows, just loop over one window
        else:
            windowlist = [Window]

        # Controls if a warning is printed in the event the background window is less than 100 data points. Hard coded to False for now
        prints_warning = False

        # Nested loops :(
        # But in reality either the middle or the outer loop is only run once
        for window in windowlist:

            # Shifts here is either the times one big background window is shifted towards the onset of the event,
            # or it is the amount of small little windows inside the larger background window
            shifts = window.n_shifts if window.n_shifts > 0 else 0

            for _ in range(shifts+1):

                # Init arrays to collect used background parameters
                mus, sigmas = np.array([]), np.array([])

                # Loop through the mus and sigmas of the current window, and find onsets 
                for j in range(window.bootstraps):

                    # Pick a number from random_choices to access a copy of flux_series in the list_of_series
                    # In the case we don't offset any of the series, there will be only one instance of flux_series in the
                    # list, and random_choices will be an array of only zeroes
                    chosen_series = list_of_series[random_choices[j]]

                    # Calculate background parameters (mean and std), and append them to their respective arrays
                    mu, sigma = sample_mean_and_std(window.start, window.end, chosen_series, sample_size=sample_size, prints_warning=prints_warning)
                    mus = np.append(mus, mu)
                    sigmas = np.append(sigmas,sigma)

                    # Choose the cusum_window to use for this iteration of the run
                    cusum_window = cusum_windows[j%len(cusum_windows)]

                    # Find an onset and save it into onset_series. Use the chosen series from the list of series
                    onset_i = onset_determination([mu, sigma], chosen_series, cusum_window, big_window_end, sigma_multiplier=sigma_multiplier)
                    onset_list.append(onset_i[-1])


                # Append the arrays to the list of all bootstrapped background params
                mus_list.append(mus)
                sigmas_list.append(sigmas)

                # Finally move the averaging window forward one hour before repeating, or in the case of small_windows,
                # add the window to the class list of windows 
                if not small_windows:
                    window.move_window(hours=1)
                else:
                    self.add_bootstrap_window(window=window)

        # After exiting the loop, return the background window to its initialized state to not cause unintended effects later
        Window.load_state(which="init")

        # Sort the list of onsets to chronological order
        onset_list.sort()

        # If de-trending is to be done, it must be done here. De-trending can only be applied if the data was resampled
        if detrend and resample:
            rs_final_index = -3 if resample[-3:] == "min" else -1
            dt_shift = f"{int(30*int(resample[:rs_final_index]))+30}s"
            onset_list = detrend_onsets(onset_list, shift=dt_shift)

        # onset_series needs to be converted to DatetimeIndex, so that the elements' mean can be calculated
        onset_indices = pd.DatetimeIndex(onset_list)
        mean_onset = onset_indices.mean()
        std_onset = np.std(onset_indices)

        # unique_vals() returns a list of the unique appearances of the onsets and their respective counts
        unique_onsets, counts = unique_vals(onset_indices)

        # First calculate the relative abundances of each unique onset
        percent_modes = np.array([count/np.nansum(counts) for count in counts])

        # The final entry of counts is the number of NaTs in the input onset_indices. Hence, it doesn't make sense
        # to keep the redundant percentage aboard. Omit the last entry to avoid errors later. No information is 
        # lost here in doing this.
        percent_modes = percent_modes[:-1]

        # Sort the percentages and unique onsets in ascending order
        percent_indices = percent_modes.argsort() # numpy array.argsort() returns a list of indices of values in ascending order

        # If absolutely no onset was found, then percent_indices will be an empty list
        if len(percent_indices) == 0:

            mean_onset = pd.NaT
            median_onset = pd.NaT
            std_onset = pd.NaT
            most_common_val = (pd.NaT, np.nan)
            onset_list = [pd.NaT]
            unique_onsets = [pd.NaT]
            confidence_intervals = [(pd.NaT, pd.NaT), (pd.NaT, pd.NaT)]

        # Here the normal case where onsets are found is treated
        else:

            unique_onsets = unique_onsets[percent_indices]
            percent_modes = percent_modes[percent_indices]

            # Unique onsets hold tuples of onset-percentage pairs in ascending order (= least likely onset first and most likely last)
            unique_onsets = np.array([(onset, percent_modes[i]) for i, onset in enumerate(unique_onsets)])

            # The most "popular" onset time and its mode are the last entries in the arrays
            most_common_val = (unique_onsets[-1,0], unique_onsets[-1,1])

            # The median onset is the middlemost onset in the distribution
            median_onset =  datetime_nanmedian(onset_list)

            # Also calculate 1-sigma and 2-sigma confidence intervals for the onset distribution
            # The method here also check the validity of the confidence intervals, which is thy it takes as an input the time resolution used
            confidence_intervals = self.get_distribution_percentiles(onset_list=onset_list, percentiles=[(15.89,84.1), (2.3,97.7)], time_reso=time_reso)

        # Attach the statistics to class attributes for easier handling. The first entry of list_of_series is the 
        # original flux_series without shifting any indices
        self.flux_series = list_of_series[0]
        bootstrap_onset_statistics = {
                                    "mean_onset" : mean_onset,
                                    "median_onset" : median_onset,
                                    "std_onset" : std_onset,
                                    "most_likely_onset" : most_common_val,
                                    "onset_list" : onset_list,
                                    "unique_onsets" : unique_onsets,
                                    "1-sigma_confidence_interval" : confidence_intervals[0],
                                    "2-sigma_confidence_interval" : confidence_intervals[1]
                                    }

        mu_and_sigma_distributions = {
                                    "mus_list" : mus_list,
                                    "sigmas_list" : sigmas_list
                                    }

        self.bootstrap_onset_statistics = bootstrap_onset_statistics
        self.mu_and_sigma_distributions = mu_and_sigma_distributions

        return bootstrap_onset_statistics, mu_and_sigma_distributions


    def show_onset_distribution(self, show_background=True, ylim=None, xlim=None, returns=False):
        """
        Displays all the unique onsets found with statistic_onset() -method in a single plot. The mode onset, that is
        the most common out of the distribution, will be shown as a solid line. All other unique onset times are shown
        in a dashed line, with the shade of their color indicating their relative frequency in the distribution.

        Note that statistic_onset() has to be run before this method. Otherwise KeyError will be raised.

        show_background : {bool}, optional
                            Boolean switch to show the background used in the plot.
        
        ylim : {list,tuple}, optional

        xlim : {list,tuple}, optional

        returns : {bool}, optional 
        """

        rcParams["font.size"] = 20

        flux_series = self.flux_series
        most_likely_onset = self.bootstrap_onset_statistics["most_likely_onset"]
        onsets = self.bootstrap_onset_statistics["unique_onsets"]

        # Create a colormap for different shades of red for the onsets. Create the map with the amount of onsets+2, so that
        # the first and the final color are left unused. They are basically white and a very very dark red.
        cmap = plt.get_cmap("Reds", len(onsets)+2)

        # Creating the figure and axes
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        # x- and y-axis settings
        ylim = set_fig_ylimits(ax=ax, flux_series=flux_series, ylim=ylim)

        ax.set_yscale("log")
        ax.set_ylabel("Intensity [1/(cm^2 sr s MeV)]")

        if not xlim:
            xlim = (flux_series.index[0], flux_series.index[-1])
        ax.set_xlim(xlim)

        set_standard_ticks(ax=ax)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(ticker.AutoLocator())
        utc_dt_format1 = DateFormatter('%H:%M \n%m-%d')
        ax.xaxis.set_major_formatter(utc_dt_format1)


        # Plots the intensity data. Keyword 'where' dictates where the step happens. Set 'mid' for the middle of the bin instead of the start.
        ax.step(flux_series.index, flux_series.values, where="mid")

        # Draws the background on the plot
        if show_background:
            self.background.draw_background(ax=ax)

        # Marking the onsets as vertical lines of varying shades of red
        for i, onset in enumerate(onsets):

            # Not the most likely onset
            if onset[0] != most_likely_onset[0]:
                linestyle = (5, (10, 3)) if onset[1] > 0.2 else (0, (5, 5))
                ax.axvline(x=onset[0], c=cmap(i+1), ls=linestyle, label=f"{str(onset[0])[11:19]} ({np.round(onset[1]*100,2):.2f} %)")

            # The most likely onset
            # I'm accessing the i+1th color value in the map because I never want the first value (that is white in the standard Reds colormap)
            else:
                ax.axvline(x=onset[0], c=cmap(i+1), label=f"{str(onset[0])[11:19]} ({np.round(onset[1]*100,2):.2f} %)") 

        # Legend and title for the figure
        ax.legend(loc=10, bbox_to_anchor=(1.0, 0.95), prop={'size': 12})
        ax.set_title("Onset distribution", fontsize=TITLE_FONTSIZE)

        plt.show()

        if returns:
            return fig, ax


    def show_onset_statistics(self, percentiles=[(15.89,84.1),(2.3,97.7)], show_background=True, xlim=None, ylim=None):
        """
        Shows the median, mode, mean and confidence intervals for the distribution of onsets got from statistic_onset()

        Parameters:
        -----------
        percentiles : {tuple, list}, optional
                    Choose which kind of percentiles of the distribution to display
        show_bacground : {bool}, optional
                    Boolean switch to show the background on the plot.
        """

        # Take the 1-std and 2-std limits from a normal distribution (default)
        sigma1, sigma2 = percentiles[0], percentiles[1]

        # Collect the 1-sigma and 2-sigma confidence intervals out of the onset distribution
        confidence_intervals = self.get_distribution_percentiles(onset_list = self.bootstrap_onset_statistics["onset_list"], percentiles=[sigma1,sigma2])

        # Collect the median, mode and mean onsets from the distribution
        onset_median = self.bootstrap_onset_statistics["median_onset"]
        onset_mode = self.bootstrap_onset_statistics["most_likely_onset"][0]
        onset_mean = self.bootstrap_onset_statistics["mean_onset"]

        # This is just for plotting with the latest resolution used
        flux = self.flux_series

        # Plot commands and settings:
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        ax.step(flux.index, flux.values, where="mid")

        # Draw the background
        if show_background:
            self.background.draw_background(ax=ax)

        # Vertical lines for the median, mode and mean of the distributions
        ax.axvline(onset_median, c="red", label="median")
        ax.axvline(onset_mode, c="blue", label="mode")
        ax.axvline(onset_mean, c="purple", label="mean")

        # 1-sigma uncertainty shading
        ax.axvspan(xmin=confidence_intervals[0][0], xmax=confidence_intervals[0][1], color="red", alpha=0.3)

        #2-sigma uncertainty shading
        ax.axvspan(xmin=confidence_intervals[1][0], xmax=confidence_intervals[1][1], color="blue", alpha=0.3)

        # Figure settings
        ax.set_yscale("log")
        ax.set_ylabel("Intensity", fontsize=AXLABEL_FONTSIZE)
        ax.set_xlabel("Time", fontsize=AXLABEL_FONTSIZE)

        ylim = set_fig_ylimits(ax=ax, flux_series=flux, ylim=ylim)

        # x-axis settings
        if not xlim:
            xlim = (flux.index[0], flux.index[-1])
        ax.set_xlim(xlim)

        set_standard_ticks(ax=ax)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(ticker.AutoLocator())
        utc_dt_format1 = DateFormatter('%H:%M \n%m-%d')
        ax.xaxis.set_major_formatter(utc_dt_format1)

        ax.legend()
        plt.show()


    def choose_flux_series(self, channels, viewing) -> tuple:
        """
        Method to choose the correct Series from a dataframe.

        Parameters:
        -----------
        channels : list
        viewing : str

        Returns:
        -----------
        flux_series : pd.Series
        en_channel_string : dict
        """

        self.choose_data(viewing)
        
        if (self.spacecraft[:2].lower() == 'st' and self.sensor == 'sept') \
                or (self.spacecraft.lower() == 'psp' and self.sensor.startswith('isois')) \
                or (self.spacecraft.lower() == 'solo' and self.sensor == 'ept') \
                or (self.spacecraft.lower() == 'solo' and self.sensor == 'het') \
                or (self.spacecraft.lower() == 'wind' and self.sensor == '3dp') \
                or (self.spacecraft.lower() == 'bepi'):
            self.viewing_used = viewing
            self.choose_data(viewing)
        elif (self.spacecraft[:2].lower() == 'st' and self.sensor == 'het'):
            self.viewing_used = ''
        elif (self.spacecraft.lower() == 'soho' and self.sensor == 'erne'):
            self.viewing_used = ''
        elif (self.spacecraft.lower() == 'soho' and self.sensor == 'ephin'):
            self.viewing_used = ''

        if self.spacecraft == 'solo':

            if self.sensor == 'het':

                if self.species in ['p', 'i']:

                    df_flux, en_channel_string =\
                        self.calc_av_en_flux_HET(self.current_df_i,
                                                 self.current_energies,
                                                 channels)
                elif self.species == 'e':

                    df_flux, en_channel_string =\
                        self.calc_av_en_flux_HET(self.current_df_e,
                                                 self.current_energies,
                                                 channels)

            elif self.sensor == 'ept':

                if self.species in ['p', 'i']:

                    df_flux, en_channel_string =\
                        self.calc_av_en_flux_EPT(self.current_df_i,
                                                 self.current_energies,
                                                 channels)
                elif self.species == 'e':

                    df_flux, en_channel_string =\
                        self.calc_av_en_flux_EPT(self.current_df_e,
                                                 self.current_energies,
                                                 channels)

            else:
                invalid_sensor_msg = "Invalid sensor!"
                raise Exception(invalid_sensor_msg)

        if self.spacecraft[:2] == 'st':

            # Super ugly implementation, but easiest to just wrap both sept and het calculators
            # in try block. KeyError is caused by an invalid channel choice.
            try:

                if self.sensor == 'het':

                    if self.species in ['p', 'i']:

                        df_flux, en_channel_string =\
                            calc_av_en_flux_ST_HET(self.current_df_i,
                                                   self.current_energies['channels_dict_df_p'],
                                                   channels,
                                                   species='p')
                    elif self.species == 'e':

                        df_flux, en_channel_string =\
                            calc_av_en_flux_ST_HET(self.current_df_e,
                                                   self.current_energies['channels_dict_df_e'],
                                                   channels,
                                                   species='e')

                elif self.sensor == 'sept':

                    if self.species in ['p', 'i']:

                        df_flux, en_channel_string =\
                            calc_av_en_flux_SEPT(self.current_df_i,
                                                 self.current_i_energies,
                                                 channels)
                    elif self.species == 'e':

                        df_flux, en_channel_string =\
                            calc_av_en_flux_SEPT(self.current_df_e,
                                                 self.current_e_energies,
                                                 channels)

            except KeyError:
                raise Exception(f"{channels} is an invalid channel or a combination of channels!")

        if self.spacecraft == 'soho':

            # A KeyError here is caused by invalid channel
            try:

                if self.sensor == 'erne':

                    if self.species in ['p', 'i']:

                        df_flux, en_channel_string =\
                            calc_av_en_flux_ERNE(self.current_df_i,
                                                 self.current_energies['channels_dict_df_p'],
                                                 channels,
                                                 species='p',
                                                 sensor='HET')

                if self.sensor == 'ephin':
                    # convert single-element "channels" list to integer
                    if type(channels) == list:
                        if len(channels) == 1:
                            channels = channels[0]
                        else:
                            print("No multi-channel support for SOHO/EPHIN included yet! Select only one single channel.")
                    if self.species == 'e':
                        df_flux = self.current_df_e[f'E{channels}']
                        en_channel_string = self.current_energies[f'E{channels}']

            except KeyError:
                raise Exception(f"{channels} is an invalid channel or a combination of channels!")

        if self.spacecraft == 'wind':
            if self.sensor == '3dp':
                # convert single-element "channels" list to integer
                if type(channels) == list:
                    if len(channels) == 1:
                        channels = channels[0]
                    else:
                        print("No multi-channel support for Wind/3DP included yet! Select only one single channel.")
                if self.species in ['p', 'i']:
                    if viewing != "omnidirectional":
                        df_flux = self.current_df_i.filter(like=f"FLUX_E{channels}_P{self.viewing[-1]}")
                    else:
                        df_flux = self.current_df_i.filter(like=f'FLUX_{channels}')
                    # extract pd.Series for further use:
                    df_flux = df_flux[df_flux.columns[0]]
                    # change flux units from '#/cm2-ster-eV-sec' to '#/cm2-ster-MeV-sec'
                    df_flux = df_flux*1e6
                    en_channel_string = self.current_i_energies['channels_dict_df']['Bins_Text'][f'ENERGY_{channels}']
                elif self.species == 'e':
                    if viewing != "omnidirectional":
                        df_flux = self.current_df_e[f"FLUX_E{channels}_P{self.viewing[-1]}"]
                        #df_flux = self.current_df_e.filter(like=f'FLUX_E{channels}')
                    else:
                        df_flux = self.current_df_e[f"FLUX_{channels}"]
                        #df_flux = self.current_df_e.filter(like=f'FLUX_{channels}')
                    # extract pd.Series for further use:
                    #df_flux = df_flux[df_flux.columns[0]]
                    # change flux units from '#/cm2-ster-eV-sec' to '#/cm2-ster-MeV-sec'
                    df_flux = df_flux*1e6
                    en_channel_string = self.current_e_energies['channels_dict_df']['Bins_Text'][f'ENERGY_{channels}']

        if self.spacecraft.lower() == 'bepi':
            if type(channels) == list:
                if len(channels) == 1:
                    # convert single-element "channels" list to integer
                    channels = channels[0]
                    if self.species == 'e':
                        df_flux = self.current_df_e[f'E{channels}']
                        en_channel_string = self.current_energies['Energy_Bin_str'][f'E{channels}']
                    if self.species in ['p', 'i']:
                        df_flux = self.current_df_i[f'P{channels}']
                        en_channel_string = self.current_energies['Energy_Bin_str'][f'P{channels}']
                else:
                    if self.species == 'e':
                        df_flux, en_channel_string = calc_av_en_flux_sixs(self.current_df_e, channels, self.species)
                    if self.species in ['p', 'i']:
                        df_flux, en_channel_string = calc_av_en_flux_sixs(self.current_df_i, channels, self.species)

        if self.spacecraft.lower() == 'psp':
            if self.sensor.lower() == 'isois-epihi':
                if self.species in ['p', 'i']:
                    # We're using here only the HET instrument of EPIHI (and not LET1 or LET2)
                    df_flux, en_channel_string =\
                        calc_av_en_flux_PSP_EPIHI(df=self.current_df_i,
                                                  energies=self.current_i_energies,
                                                  en_channel=channels,
                                                  species='p',
                                                  instrument='het',
                                                  viewing=viewing.upper())
                if self.species == 'e':
                    # We're using here only the HET instrument of EPIHI (and not LET1 or LET2)
                    df_flux, en_channel_string =\
                        calc_av_en_flux_PSP_EPIHI(df=self.current_df_e,
                                                  energies=self.current_e_energies,
                                                  en_channel=channels,
                                                  species='e',
                                                  instrument='het',
                                                  viewing=viewing.upper())
            if self.sensor.lower() == 'isois-epilo':
                if self.species == 'e':
                    # We're using here only the F channel of EPILO (and not E or G)
                    df_flux, en_channel_string =\
                        calc_av_en_flux_PSP_EPILO(df=self.current_df_e,
                                                  en_dict=self.current_e_energies,
                                                  en_channel=channels,
                                                  species='e',
                                                  mode='pe',
                                                  chan='F',
                                                  viewing=viewing)

        if self.spacecraft == 'solo':
            flux_series = df_flux #[channels]
        if self.spacecraft[:2].lower() == 'st':
            flux_series = df_flux  # [channel]'
        if self.spacecraft.lower() == 'soho':
            flux_series = df_flux  # [channel]
        if self.spacecraft.lower() == 'wind':
            flux_series = df_flux  # [channel]
        if self.spacecraft.lower() == 'psp':
            flux_series = df_flux #[channels]
        if self.spacecraft.lower() == 'bepi':
            flux_series = df_flux  # [channel]


        # Before returning, make sure that the type is pandas series, and not a 1-dimensional dataframe
        if not isinstance(flux_series, pd.core.series.Series):
            flux_series = flux_series.squeeze()

        return flux_series, en_channel_string


    def plot_subwindows(self):
        """
        Shows the subwindows used in the bootstrapping process
        """

        try:
            subwindows_plotted = self.subwindows_plotted
        except AttributeError:
            print("The onset object currently has no subwindows!")
            return 0

        fig, ax = self.fig, self.ax

        if not self.subwindows_plotted:
            for i, window in enumerate(self.list_of_bootstrap_windows):

                # Color the window
                ax.axvspan(xmin=window.start, xmax=window.end,
                            color=self.window_colors[i], alpha=0.25)

                # Draw the window boundaries with thin black lines
                ax.axvline(x=window.start, c='k', lw=1)
                ax.axvline(x=window.end, c='k', lw=1)

            # Keep track of wether subwindows have been plotted, to avoid over-plotting
            self.subwindows_plotted = True

        #plt.show()
        return fig


    def plot_cdf(self, prints=False) -> None:
        """
        Plots the Cumulatice Distribution Function (CDF) for a distribution of onset times.

        Parameters:
        -----------
        prints : bool, default False
                If True, will also print out the median onset time and the confidence intervals.
        """

        # Remember that onset times are sorted according to their probability in ascending order
        onset_timestamps = np.array([pair[0] for pair in self.bootstrap_onset_statistics["unique_onsets"]])
        probabilities = np.array([pair[1] for pair in self.bootstrap_onset_statistics["unique_onsets"]])

        # Here let's sort the onsets and probabilities in temporal order:
        onset_time_indices = onset_timestamps.argsort() # numpy array.argsort() returns a list of indices of values in ascending order

        onset_timestamps = pd.to_datetime(onset_timestamps[onset_time_indices])
        probabilities = probabilities[onset_time_indices]

        # Effectively calculate the cumulative sum of probabilities
        cumulative_probabilities = [np.nansum(probabilities[:i+1]) for i in range(len(probabilities))]
        
        # Get the median, mode, mean and confidence intervals (1sigma and 2sigma of normal distribution) of the onset distribution
        median_onset, _, _, confidence_intervals_all = self.get_distribution_statistics(onset_statistics=self.bootstrap_onset_statistics,
                                                                                    percentiles=[(15.89,84.1), (2.3,97.7)])

        # Flatten the sigma lists for plotting. Look at this list comprehension :D
        confidence_intervals_all = [timestamp for sublist in confidence_intervals_all for timestamp in sublist]

        # Init the figure
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        rcParams["font.size"] = 20

        # Settings for axes
        ax.xaxis_date()
        ax.set_xlabel("Time", fontsize=AXLABEL_FONTSIZE)
        ax.set_ylabel("Cumulative Probability", fontsize=AXLABEL_FONTSIZE)
        ax.set_title("Cumulative Distribution Function", fontsize=TITLE_FONTSIZE)
        hour_minute_format = DateFormatter("%H:%M")
        ax.xaxis.set_major_formatter(hour_minute_format)

        # Plot command(s)
        ax.step(onset_timestamps, cumulative_probabilities, zorder=3, c='k', where="post")

        colors = ("red", "red", "blue", "blue")
        # Horizontal dashed lines for 1-sigma and 2-sigma limits
        for i, num in enumerate((0.159, 0.841, 0.023, 0.977)):
            label = "15.9-84.1 percentiles" if num==0.841 else ("2.3-97.7 percentiles" if num==0.977 else None)
            ax.axhline(y=num, ls="--", c=colors[i], zorder=1, alpha=0.8, label=label)

        # Vertical dashed lines for 1-sigma and 2-sigma limits
        for i, time in enumerate(confidence_intervals_all):
            ax.axvline(x=time, ls="--", c=colors[i], zorder=1, alpha=0.8)

        ax.axvline(x=median_onset, ls="--", c="maroon", zorder=2, alpha=0.8, label="median onset (50th percentile)")
        ax.axhline(y=0.50, ls="--", c="maroon", zorder=1, alpha=0.8)
        
        ax.legend(loc=10, bbox_to_anchor=(1.27, 0.9))

        if prints:
            print(f"The median onset time is: {str(median_onset.time())[:8]}")
            # print(f"The np.percentile()-calculated median onset time (50th percentile) is: {str(median_onset.time())[:8]}")
            display(Markdown("~68 % confidence interval"))
            print(f"{str(confidence_intervals_all[0].time())[:8]} - {str(confidence_intervals_all[1].time())[:8]}")
            print(" ")
            display(Markdown("~95 % confidence interval"))
            print(f"{str(confidence_intervals_all[2].time())[:8]} - {str(confidence_intervals_all[3].time())[:8]}")

        plt.show()


    def scatter_histogram(self, x:str="mean", xbinwidth:int=1, ybinwidth:str="1min") -> None:
        """
        A method to plot a scatter and histogram plots of either background mean or background std
        vs. onset time.

        Parameters:
        -----------
        x : str, default 'mean'
            Either 'mean' or 'std'
        xbinwidth : int
            The width of x-axis bins
        ybinwidth: str
            Pandas-compatible time string for the width of onset distribution.
        """

        x_axes_choice = {
            "mean" : self.mu_and_sigma_distributions["mus_list"],
            "std" : self.mu_and_sigma_distributions["sigmas_list"]
        }

        xdata = x_axes_choice[x]
        ydata = self.bootstrap_onset_statistics["onset_list"]

        title = r"Background $\mu$ vs onset time" if x=="mean" else r"Background $\sigma$ vs onset time"
        xlabel = r"Sample mean intensity" if x=="mean" else r"Sample standard deviation"
        ylabel = "Onset time [D HH:MM]"

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        spacing = 0.005

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]

        rcParams["font.size"] = 20

        def scatter_hist(x, y, ax, ax_histx, ax_histy) -> None:

            # No labels for the histograms
            ax_histx.tick_params(axis="x", labelbottom=False)
            ax_histy.tick_params(axis="y", labelleft=False)

            # The scatter plot:
            for i, arr in enumerate(x):
                element = len(arr)
                ax.scatter(arr, y[i*element:element*(1+i)], color=self.window_colors[i], s=10)

            ylims = (np.min(y) - pd.Timedelta(minutes=2), np.max(y) + pd.Timedelta(minutes=2))
            ax.set_ylim(ylims)
            ax.grid(True)

            # Now determine limits by hand:
            binwidth = xbinwidth
            xymax = int(np.max(x))
            xmin = int(np.min(x))

            xbins = np.arange(xmin, xymax, binwidth)
            for i, arr in enumerate(x):
                ax_histx.hist(arr, bins=xbins, color=self.window_colors[i], alpha=0.6)

            half_bin = pd.Timedelta(seconds=30)
            ybins = pd.date_range(start=ylims[0]+half_bin, end=ylims[1]+half_bin, freq=ybinwidth).tolist()

            onset_frequencies = np.ones_like(y)/len(y)

            ax_histy.hist(y, bins=ybins, edgecolor="black", orientation='horizontal', weights=onset_frequencies)

            max_freq = np.nanmax([pair[1] for pair in self.bootstrap_onset_statistics["unique_onsets"]])
            ax_histy.set_xticks([np.round(max_freq,2)/4, np.round(max_freq,2)/2, np.round(max_freq,2)])


        # Start with a square-shaped Figure
        fig = plt.figure(figsize=(12, 12))

        # These axes are for the scatter plot in the middle
        ax = fig.add_axes(rect_scatter)
        ax.yaxis_date()
        
        # These axes are for the histograms
        ax_histx = fig.add_axes(rect_histx, sharex=ax)
        ax_histy = fig.add_axes(rect_histy, sharey=ax)

        # Use the function defined above
        scatter_hist(xdata, ydata, ax, ax_histx, ax_histy)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        fig.suptitle(title, fontsize=20)

        # plt.savefig("histogram_bg_vars")
        plt.show()


    def VDA(self, onset_times=None, Onset=None, energy:str='gmean', selection=None, 
            yerrs=None, reference:str="mode", title=None, ylim=None, plot=True, guess=None, save=False,
            savepath=None, grid=True, show_omitted=True):
        """
        Performs Velocity Dispersion Analysis.

        Parameters:
        -----------
        onset_times : array-like, optional
                    List of onset times. If None, automatically acquired from self.onset_statistics
        Onset : self@Onset
                    The second onset object if joint-VDA
        energy: str, default 'gmean'
                    Mode of assessing nominal channel energy. 'gmean' or 'upper'
        selection: slice(start,stop) or a pair of slices, default: None
                    Selection of the channels one wishes to include in VDA calculation. If a pair of slices are delivered,
                    consider each slice to apply individually to each instrument respectively.
        yerrs : {array-like}, optional
                    Use provided y-errors instead of the ones in the object's memory.
        reference : {str}, either 'median' or 'mode' (default)
                    Which one to consider as the onset times extracted with the statistical method.
        title : {str}, default None
                    Custom title. If None, generate standard title.
        ylim : tuple/list with len==2, default None
                    the limits of y-axis.
        plot: {bool}, default True
                    Switch to produce plot of the analysis.
        guess : {array-like} of len()==2, default None
                    The initial guess values for the ODR fit. First one's the slope, which is a float, and the second is the
                    intersection point, which is a pandas-compatible datetime string.
        save: {bool}, default False
                    Switch to save the plotted figure. Only works if plot=True.
        grid : {bool} default True
                    Boolean switch for gridlines.
        show_omitted : {bool} default True
                    A switch to show the data points omitted from the fit in the plot.

        Returns:
        ---------
        output: dict
                    A dictionary containing  inverse_betas, x_errors, onset_times, y_errors, path_length, injection_time,
                    path_length uncertainty, injection_time uncertainty, residual variance (ODR),
                    stopreason (ODR), fig and axes.
        """

        VDA_MARKERSIZE = 24
        VDA_ELINEWIDTH = 2.2
        VDA_CAPSIZE = 6.0

        from scipy.stats import t as studentt

        import numpy.ma as ma

        spacecraft = self.spacecraft.lower()
        instrument = self.sensor.lower()
        species = self.species.lower()

        if species in ELECTRON_IDENTIFIERS:
            species_title = "electrons"
            m_species = const.m_e.value
        elif species in ("ion", 'i', 'h'):
            species_title = "ions"
            m_species = const.m_p.value
        elif species in ("protons", "proton", 'p'):
            species_title = "protons"
            m_species = const.m_p.value
        else:
            raise Exception(f"Particle species '{self.species}' does not appear to be any of the recognized species!")

        # E=mc^2, a fundamental property of an object with mass
        mass_energy = m_species*C_SQUARED # ~511 keV for electrons

        # Get the energies of each energy channel, to calculate the mean energy of particles and ultimately
        # to get the dimensionless speeds of the particles (v/c = beta).
        # This method returns lower and higher energy bounds in electron volts
        if self.spacecraft in SEPPY_SPACECRAFT:
            e_min, e_max = self.get_channel_energy_values()
        else:
            e_min, e_max = self.get_custom_channel_energies()

         # Initialize the list of channels according to the sc and instrument used
        channels = self.get_all_channels()

        # SOHO /EPHIN really has only 4 active channels, but 5th one still exists, causing an erroneous amount
        # of channel nominal energies and hence wrong sized mask (5 instead of 4). For now we just discard
        # the final entries fo these lists.
        if spacecraft=="soho" and instrument=="ephin" and len(e_min)==5:
            e_min, e_max = e_min[:-1], e_max[:-1]

        # Get the nominal channel energies, which are by default the geometric means of channel boundaries.
        nominal_energies = calc_chnl_nominal_energy(e_min, e_max, mode=energy)

        # Check here if onset_times were given as an input. If not, use median/mode onset times and yerrs.
        onset_times = np.array([]) if onset_times is None else np.array(onset_times)
        if len(onset_times)==0:

            # Include a check here to get rid of channel 300 in EPHIN data after the switch-off
            if self.sensor == "ephin" and self.start_date > EPHIN_300_INVALID_ONWARDS:
                channels = [c_id for c_id in channels if c_id != 300]

            for ch in channels:
                # The median onset is found as the first entry of any list within the dictionary 'onset_statistics'. If that,
                # however, does not exist for a particular channel, then insert NaT to be masked away later in the method.
                try:
                    onset_times = np.append(onset_times,self.onset_statistics[ch][0])
                except KeyError as e:
                    onset_times = np.append(onset_times,pd.NaT)

        # Here check if two lists of onset times are being performed VDA upon
        if Onset:

            # The other set of onset times may belomng to a different particle species
            if Onset.species in ELECTRON_IDENTIFIERS:
                m_species1, species_title1 = const.m_e.value, "electrons"
            else:
                m_species1, species_title1 = const.m_p.value, "protons"
            mass_energy1 = m_species1 * C_SQUARED

            # First get the energies and nominal energies for the second instrument
            if Onset.spacecraft in SEPPY_SPACECRAFT:
                e_min1, e_max1 = Onset.get_channel_energy_values()
            else:
                e_min1, e_max1 = Onset.get_custom_channel_energies()

            # SOHO /EPHIN really has only 4 active channels, but 5th one still exists, causing an erroneous amount
            # of channel nominal energies and hence wrong sized mask (5 instead of 4). For now we just discard
            # the final entries fo these lists.
            if Onset.spacecraft=="soho" and Onset.sensor=="ephin" and len(e_min1)==5:
                e_min1, e_max1 = e_min1[:-1], e_max1[:-1]

            nominal_energies1 = calc_chnl_nominal_energy(e_min1, e_max1, mode=energy)

             # Check here if onset_times were given as an input. If not, use median/mode onset times and yerrs.
            if len(onset_times) == 2 and isinstance(onset_times[0], (list,tuple, np.ndarray)):

                onset_times1 = pd.to_datetime(onset_times[1])
                onset_times = pd.to_datetime(onset_times[0])

            else:

                # initialize the list of channels according to the sc and instrument used
                channels1 = Onset.get_all_channels()

                # Include a check here to get rid of channel 300 in EPHIN data after the switch-off
                if Onset.sensor == "ephin" and Onset.start_date > EPHIN_300_INVALID_ONWARDS:
                    channels1 = [c_id for c_id in channels1 if c_id != 300]

                # If there was no input for onset times, init an empty array and fill it up with values from the object
                onset_times1 = np.array([])

                for ch in channels1:

                    # The mode onset is found as the first entry of any list within the dictionary 'onset_statistics'. If that,
                    # however, does not exist for a particular channel, then insert NaT to be masked away later in the method.
                    try:
                        onset_times1 = np.append(onset_times1,Onset.onset_statistics[ch][0])
                    except KeyError as e:
                        onset_times1 = np.append(onset_times1,pd.NaT)


            onset_times_all = np.concatenate((onset_times, onset_times1))


            # Calculate the inverse betas corresponding to the nominal channel energies
            inverse_beta = calculate_inverse_betas(channel_energies=nominal_energies, mass_energy=mass_energy)
            inverse_beta1 = calculate_inverse_betas(channel_energies=nominal_energies1, mass_energy=mass_energy1)

            inverse_beta_all = np.concatenate((inverse_beta, inverse_beta1))

            # This is for the fitting function. 8.33 min = light travel time/au
            inverse_beta_corrected = np.array([8.33*60*b for b in inverse_beta_all]) #multiply by 60 -> minutes to seconds

            # Second values of the onset times
            date_in_sec = datetime_to_sec(onset_times=onset_times)
            date_in_sec1 = datetime_to_sec(onset_times=onset_times1)

            date_in_sec_all = np.concatenate((date_in_sec, date_in_sec1))

            # Error bars in x direction:
            x_errors_lower, x_errors_upper,  x_errors = get_x_errors(e_min=e_min, e_max=e_max, inverse_betas=inverse_beta, mass_energy=mass_energy)
            x_errors_lower1, x_errors_upper1,  x_errors1 = get_x_errors(e_min=e_min1, e_max=e_max1, inverse_betas=inverse_beta1, mass_energy=mass_energy1)

            # Arrays to hold all the lower and upper energy bounds. These might be redundant
            x_errors_lower_all = np.concatenate((x_errors_lower, x_errors_lower1))
            x_errors_upper_all = np.concatenate((x_errors_upper, x_errors_upper1))

            x_errors_all = np.concatenate((x_errors, x_errors1))

            # Loop through all possible, channels, even those that not necessarily show an onset
            # Remember:  self.onset_statistics : {channel_id : [mode, median, 1st_sigma_minus, 1st_sigma_plus, 2nd_sigma_minus, 2nd_sigma_plus]}
            if reference == "mode":
                ref_idx = 0
            elif reference == "median":
                ref_idx = 1
            else:
                raise ValueError(f"Argument {reference} is an invalid input for the variable 'reference'. Acceptable input values are are 'mode' and 'median'.")

            # Get all the y-errors there are in this object's database
            if not isinstance(yerrs, (list,np.ndarray)):

                plus_errs, minus_errs = np.array([]), np.array([])
                plus_errs1, minus_errs1 = np.array([]), np.array([])

                for ch in channels:

                    try: 
                        minus_err = self.onset_statistics[ch][ref_idx] - self.onset_statistics[ch][4] # the difference of a timestamp and a timedelta is a timedelta
                        plus_err = self.onset_statistics[ch][5] - self.onset_statistics[ch][ref_idx] 
                    except KeyError as e:
                        print(f"KeyError in channel {e}. Missing onset?")
                        # It's irrelevant what timedelta we insert here; the corresponding data point does not exist
                        plus_err = pd.Timedelta(seconds=1)
                        minus_err = pd.Timedelta(seconds=1)

                    # Before a check was done here to make sure that the errors are not smaller than half of the minimum
                    # cadence of the respective instrument; this is now obsolete, as the errors are properly inspected
                    # while the weighting is performed.
                    # plus_errs = np.append(plus_errs, plus_err) if plus_err >= self.minimum_cadences[f"{spacecraft}_{instrument}"]/2 else np.append(plus_errs, pd.Timedelta(self.minimum_cadences[f"{spacecraft}_{instrument}"])/2)
                    # minus_errs = np.append(minus_errs, minus_err) if minus_err >= self.minimum_cadences[f"{spacecraft}_{instrument}"]/2 else np.append(minus_errs, pd.Timedelta(self.minimum_cadences[f"{spacecraft}_{instrument}"])/2)
                    plus_errs = np.append(plus_errs, plus_err)
                    minus_errs = np.append(minus_errs, minus_err)


                for ch in channels1:

                    try: 
                        minus_err1 = Onset.onset_statistics[ch][ref_idx] - Onset.onset_statistics[ch][4] # the difference of a timestamp and a timedelta is a timedelta
                        plus_err1 = Onset.onset_statistics[ch][5] - Onset.onset_statistics[ch][ref_idx] 
                    except KeyError as e:
                        print(f"KeyError in channel {e}. Missing onset?")
                        # It's irrelevant what timedelta we insert here; the corresponding data point does not exist
                        plus_err1 = pd.Timedelta(seconds=1)
                        minus_err1 = pd.Timedelta(seconds=1)

                    plus_errs1 = np.append(plus_errs1, plus_err1)
                    minus_errs1 = np.append(minus_errs1, minus_err1)

                plus_errs_all = np.append(plus_errs, plus_errs1)
                minus_errs_all = np.append(minus_errs, minus_errs1)

                # Convert errors in time to seconds
                plus_errs_secs = [err.seconds for err in plus_errs]
                minus_errs_secs = [err.seconds for err in minus_errs]
                plus_errs_secs1 = [err.seconds for err in plus_errs1]
                minus_errs_secs1 = [err.seconds for err in minus_errs1]

                plus_errs_secs_all = np.append(plus_errs_secs, plus_errs_secs1)
                minus_errs_secs_all = np.append(minus_errs_secs, minus_errs_secs1)

                # Uneven errorbars need to be shape (2,N), where first row contains the lower errors, the second row contains the upper errors.
                y_errors_all = np.array([minus_errs_secs_all, plus_errs_secs_all])

            else:
                # The y-directional error was given by the user
                # 4 arrays, so asymmetric plus and minus errors
                if len(yerrs)==4:

                    plus_errs, minus_errs = np.array([]), np.array([])
                    plus_errs1, minus_errs1 = np.array([]), np.array([])

                    minus_timestamps, plus_timestamps = yerrs[0], yerrs[1]
                    minus_timestamps1, plus_timestamps1 = yerrs[0], yerrs[1]

                    for i, ch in enumerate(channels):

                        try: 
                            minus_err = minus_timestamps[i]
                            plus_err = plus_timestamps[i]
                        except KeyError as e:
                            plus_err = pd.Timedelta(seconds=1)
                            minus_err = pd.Timedelta(seconds=1)
                        
                        plus_errs = np.append(plus_errs, plus_err)
                        minus_errs = np.append(minus_errs, minus_err)

                    y_errors_all = np.array([plus_errs,minus_errs])

                    # Convert errors in time to seconds
                    plus_errs_secs = [err.seconds for err in plus_errs]
                    minus_errs_secs = [err.seconds for err in minus_errs]

                    y_errors_all_secs = np.array([plus_errs_secs,minus_errs_secs])

                # 2 arrays -> symmetric errors for both
                elif len(yerrs)==2:

                    # Check the first entry of yerrs. If it's not a timedelta, then it's the reach of the error
                    if not isinstance(yerrs[0][0], (pd.Timedelta, datetime.timedelta)):
                        y_errors_all = np.array([])  #= yerrs

                        for i, ch in enumerate(channels):
                            y_err = yerrs[i] - self.onset_statistics[ch][ref_idx]
                            y_errors_all = np.append(y_errors, y_err)

                    else:
                        # These are for the fitting function
                        y_errors = yerrs[0]
                        y_errors1 = yerrs[1]

                        y_errors_all = np.concatenate((y_errors, y_errors1))
                        y_errors_all_secs = [err.seconds for err in y_errors_all]

                # 1 array -> wrong
                else:
                    raise ValueError("The y-errors for two-instrument VDA must be in form of 4 or 2 lists!")

            # Numpy masks work so that True values get masked as invalid, while False remains unaltered
            mask = np.isnan(date_in_sec_all)

            # These can be used as maskedarray even when no invalid values are in the onset times
            date_in_sec_all = ma.array(date_in_sec_all, mask=mask)
            onset_times_all = ma.array(onset_times_all, mask=mask)
            inverse_beta_all = ma.array(inverse_beta_all, mask=mask)
            inverse_beta_corrected = ma.array(inverse_beta_corrected, mask=mask)
            x_errors_all = ma.array(x_errors_all, mask=mask)

            # Asymmetric errors / errors not defined by the user
            if len(y_errors_all)==2:
                y_errors_plot = ma.array([minus_errs, plus_errs], mask=[np.isnan(date_in_sec),np.isnan(date_in_sec)])
                y_errors1_plot = ma.array([minus_errs1, plus_errs1], mask=[np.isnan(date_in_sec1),np.isnan(date_in_sec1)])
                y_errors_all_plot = ma.array([minus_errs_all, plus_errs_all], mask=[mask,mask])
                y_errors_all_secs = ma.array([minus_errs_secs_all, plus_errs_secs_all], mask=[mask,mask])

            else:
                y_errors_plot = ma.array(y_errors, mask=np.isnan(date_in_sec))
                y_errors1_plot = ma.array(y_errors1, mask=np.isnan(date_in_sec1))
                y_errors_all_plot = ma.array(y_errors_all, mask=mask)
                y_errors_all_secs = ma.array(y_errors_all_secs, mask=mask)

            # After masking NaNs and Nats away, slice
            # which datapoints to consider for the fit
            if selection is not None and isinstance(selection[0],(slice,list)):

                if isinstance(selection[0],slice):
                    selection1 = selection[1]
                    selection = selection[0]

                    # The joint selection is a little tricky, as it requires concatenating two selections
                    bool_mask = [True if (i >= selection.start and i < selection.stop) else False for i in range(len(onset_times))]

                    # Check selection1 type to treat accordingly
                    if isinstance(selection1, slice):
                        bool_mask1 = [True if (i >= selection1.start and i < selection1.stop) else False for i in range(len(onset_times1))]
                    else:
                        bool_mask1 = selection1

                    # Finally combine selections
                    selection_all = np.append(bool_mask, bool_mask1)

                # In this case selection is a boolean [TRUE TRUE FALSE ...] list
                else:

                    selection1 = selection[1]
                    selection = selection[0]

                    # Check selection1 type to treat accordingly
                    if isinstance(selection1, slice):
                        selection1 = [True if (i >= selection1.start and i < selection1.stop) else False for i in range(len(onset_times1))]

                    # Finally combine selections
                    selection_all = np.append(selection,selection1) 

                # This variable signifies if there are data points in the plot that are not considered for the fit,
                # and are therefore hollowed with a white middle part.
                omitted_exists = True

            else:
                selection = slice(0,len(onset_times))
                selection1 = slice(0,len(onset_times1))
                selection_all = slice(0,len(onset_times_all))
                omitted_exists = False

            # These are only used for the fit -> slice them to fit selection
            inverse_beta_corrected = inverse_beta_corrected[selection_all]
            date_in_sec_all = date_in_sec_all[selection_all]

        # Only one instrument:
        else:

            # Calculate the inverse betas corresponding to the nominal channel energies
            inverse_beta = calculate_inverse_betas(channel_energies=nominal_energies, mass_energy=mass_energy)

            # This is for the fitting function. 8.33 min = light travel time/au , coeff = 8.33*60
            inverse_beta_corrected = np.array([b*8.33*60 for b in inverse_beta]) #multiply by 60 -> minutes to seconds

            # Get the second values for onset times for the fit
            date_in_sec = datetime_to_sec(onset_times=onset_times)

            # Error bars in x direction:
            x_errors_lower, x_errors_upper,  x_errors_all = get_x_errors(e_min=e_min, e_max=e_max, inverse_betas=inverse_beta, mass_energy=mass_energy)


            # Choose the reference index here before calculating the timedeltas for y-errors
            # self.onset_statistics : {channel_id : [mode, median, 1st_sigma_minus, 1st_sigma_plus, 2nd_sigma_minus, 2nd_sigma_plus]}
            if reference == "mode":
                ref_idx = 0
            elif reference == "median":
                ref_idx = 1
            else:
                raise ValueError(f"Argument {reference} is an invalid input for the variable 'reference'. Acceptable input values are are 'mode' and 'median'.")

            # Get all the y-errors there are in this object's database
            if not isinstance(yerrs, (list, np.ndarray)):

                plus_errs, minus_errs = np.array([]), np.array([])
                # Loop through all possible, channels, even those that not necessarily show an onset
                for ch in channels:

                    try: 
                        minus_err = self.onset_statistics[ch][ref_idx] - self.onset_statistics[ch][4] # the difference of a timestamp and a timedelta is a timedelta
                        plus_err = self.onset_statistics[ch][5] - self.onset_statistics[ch][ref_idx] 
                    except KeyError as e:
                        plus_err = pd.Timedelta(seconds=1)
                        minus_err = pd.Timedelta(seconds=1)

                    # The previous way of determining the validity of the ~95% confidence intervals here should be redundant now that the widths of 
                    # the intervals and their boundaries are already checked when they are determined. Instead, just collect the errors to the arrays.
                    # plus_errs = np.append(plus_errs, plus_err) if plus_err >= self.minimum_cadences[f"{spacecraft}_{instrument}"]/2 else np.append(plus_errs, pd.Timedelta(self.minimum_cadences[f"{spacecraft}_{instrument}"])/2)
                    # minus_errs = np.append(minus_errs, minus_err) if minus_err >= self.minimum_cadences[f"{spacecraft}_{instrument}"]/2 else np.append(minus_errs, pd.Timedelta(self.minimum_cadences[f"{spacecraft}_{instrument}"])/2)
                    plus_errs = np.append(plus_errs, plus_err)
                    minus_errs = np.append(minus_errs, minus_err)

                # Convert errors in time to seconds
                plus_errs_secs = [err.seconds for err in plus_errs]
                minus_errs_secs = [err.seconds for err in minus_errs]

                # Uneven errorbars need to be shape (2,N), where first row contains the lower errors, the second row contains the upper errors.
                y_errors_all = np.array([minus_errs, plus_errs])

            # User gave y-errors as an input
            else:

                # 2 arrays, so asymmetric plus and minus errors
                if len(yerrs)==2:

                    plus_errs, minus_errs = np.array([]), np.array([])

                    minus_timestamps, plus_timestamps = yerrs[0], yerrs[1]

                    for i, ch in enumerate(channels):

                        try: 
                            minus_err = minus_timestamps[i]
                            plus_err = plus_timestamps[i]
                        except KeyError as e:
                            plus_err = pd.Timedelta(seconds=1)
                            minus_err = pd.Timedelta(seconds=1)
                        
                        plus_errs = np.append(plus_errs, plus_err)
                        minus_errs = np.append(minus_errs, minus_err)

                    y_errors_all = np.array([plus_errs,minus_errs])

                    # Convert errors in time to seconds
                    plus_errs_secs = [err.seconds for err in plus_errs]
                    minus_errs_secs = [err.seconds for err in minus_errs]

                    y_errors_all_secs = np.array([plus_errs_secs,minus_errs_secs])

                # 1 array -> symmetric errors
                else:

                    # Check the first entry of yerrs. If it's not a timedelta, then it's the reach of the error
                    if not isinstance(yerrs[0], (pd.Timedelta, datetime.timedelta)):
                        y_errors_all = np.array([])  #= yerrs

                        for i, ch in channels:
                            y_err = yerrs[i] - self.onset_statistics[ch][ref_idx]
                            y_errors_all = np.append(y_errors, y_err)

                    else:
                        y_errors_all = yerrs

                    y_errors_all_secs = [err.seconds for err in y_errors_all]


            # Numpy masks work so that True values get masked as invalid, while False remains unaltered
            mask = np.isnan(date_in_sec)

            # These can be used as maskedarray even when no invalid values are in the onset times
            date_in_sec = ma.array(date_in_sec, mask=mask)
            onset_times = ma.array(onset_times, mask=mask)
            inverse_beta = ma.array(inverse_beta, mask=mask)
            inverse_beta_corrected = ma.array(inverse_beta_corrected, mask=mask)
            x_errors_all = ma.array(x_errors_all, mask=mask)

            # Errors are asymmetric
            if len(y_errors_all)==2:
                y_errors_plot = ma.array([minus_errs, plus_errs], mask=[mask,mask])
                y_errors_all_secs = ma.array([minus_errs_secs, plus_errs_secs], mask=[mask,mask])

            # Errors are symmetric
            else:
                y_errors_all_secs = ma.array(y_errors_all_secs, mask=mask)
                y_errors_plot = ma.array(y_errors_all, mask=mask)

            # Just to have a shared alias for all the y-errors, also in the case there are two instruments
            y_errors_all_plot = y_errors_plot

            # After masking NaNs and Nats away, slice
            # which datapoints to consider for the fit
            if selection is not None and isinstance(selection,(slice,list)):

                # User may give empty slice, e.g., slice(0,0) as an input. Handle it as if selection was None.
                if isinstance(selection,slice):
                    if selection.start==selection.stop:
                        selection_all = ~mask
                        selection = selection_all
                        omitted_exists = False
                
                selection_all = selection
                omitted_exists = True

            else:
                selection_all = ~mask
                selection = selection_all
                omitted_exists = False
                

            # These are only used for the fit -> slice them to fit selection
            inverse_beta_corrected = inverse_beta_corrected[selection_all]
            date_in_sec_all = date_in_sec[selection_all]

            # Common name to take into account single instrument and two-instrument code blocks
            inverse_beta_all = inverse_beta
            onset_times_all = onset_times

        # Here onward we do the fit, and find the slope and intersection of this fit ->

        # Here happens the fit and errors ->
        # x_errors and y_errors are fed in here with selections because we want to plot them for all data points, but only
        # consider the selection of them in the fit
        # The .compressed() -method here ensures that only the valid values of the masked arrays are fed into the function. If this is not
        # done, then all of the values that the function yields will be nan.
        odr_output = seek_fit_and_errors(x=inverse_beta_corrected.compressed(), y=date_in_sec_all.compressed(), 
                                         xerr=x_errors_all.compressed(), yerr=y_errors_all_secs.compressed(), guess=guess)
        
        # This would print out the output
        # odr_output.pprint()
        
        slope, constant = odr_output.beta[0], odr_output.beta[1]
        
        dof = len(inverse_beta_corrected) - len(odr_output.beta) # degrees of freedom (number of data points minus number of free parameters)
        t_val = studentt.interval(0.95, dof)[1] # 95 % confidence interval
        errors = t_val * odr_output.sd_beta
        slope_uncertainty = errors[0]

        try:
            t_inj_uncertainty = pd.Timedelta(seconds=int(errors[1]))

        # Valueerror is caused by nan
        except ValueError:
            t_inj_uncertainty = np.nan

        # Residual variance
        residual_variance = odr_output.res_var

        # Reason for stopping the regression:
        stopreason = odr_output.stopreason

        odr_fit_sec = inverse_beta_corrected*slope
        odr_fit = np.array([datetime.datetime.utcfromtimestamp(sec + constant) for sec in odr_fit_sec.compressed()])

        # Release time is the timestamp where the line intersects y-axis
        release_time = datetime.datetime.utcfromtimestamp(constant)

        # Precision of t_inj display:
        precision = 5 if spacecraft not in ("solo","wind") else 8
        rel_time_str =str(release_time.time())[:precision]

        date_of_event = get_figdate(onset_times_all.compressed())

        # Declare axes, the things plotted and returned by this function
        # First observed datapoints, then fit and its slope+constant, lastly y-and x errors

        output = {
                "inverse_beta" : inverse_beta_all,
                "x_errors" : x_errors_all,
                "onset_times" : onset_times_all,
                "y_errors" : y_errors_all_plot,
                "path_length" : slope,
                "t_inj" : release_time,
                "path_length_uncertainty" : slope_uncertainty,
                "t_inj_uncertainty" : t_inj_uncertainty,
                "residual_variance" : residual_variance,  # The unexplained variance of the data according to the model
                "stopreason" : stopreason
                 }

        # Only plotting commands from here onward->
        if plot:

            plt.rcParams['axes.linewidth'] = 2.5
            plt.rcParams['font.size'] = 16

            # The colorblind style is apparently dependent on the configuration of stars and planets,
            # because sometimes 'seaborn-v0_8-colorblind works and sometimes it doesn't. So let's
            # try/except here to find a style that actually works.
            try:
                plt.style.use("seaborn-v0_8-colorblind")
            except OSError:
                plt.style.use("seaborn-colorblind")

            fig, ax = plt.subplots(figsize=VDA_FIGSIZE)

            ax.grid(visible=grid, axis="both")

            # Choose if to plot all data points, or just the ones that are considered for the fit. The default is to show all data points.
            if not show_omitted:
                omitted_exists = False
                inverse_beta = inverse_beta[selection]
                onset_times = onset_times[selection]
                x_errors_lower = x_errors_lower[selection]
                x_errors_upper = x_errors_upper[selection]
                y_errors_plot = y_errors_plot.compressed()
                if len(y_errors_plot)==2:
                    y_errors_plot[0] = y_errors_plot[0][selection]
                    y_errors_plot[1] = y_errors_plot[1][selection]
                else:
                    y_errors_plot = y_errors_plot[selection]

                if Onset:
                    inverse_beta1 = inverse_beta1[selection1]
                    onset_times1 = onset_times1[selection1]
                    x_errors_lower1 = x_errors_lower1[selection1]
                    x_errors_upper1 = x_errors_upper1[selection1]
                    if len(y_errors1_plot)==2:
                        y_errors1_plot_tmp0 = y_errors1_plot[0][selection1]
                        y_errors1_plot_tmp1 = y_errors1_plot[1][selection1]
                        y_errors1_plot = [y_errors1_plot_tmp0, y_errors1_plot_tmp1]
                    else:
                        y_errors1_plot = y_errors1_plot[selection1]


            # About matplotlib.Axes.errorbar:
            # shape(2, N): Separate - and + values for each bar. First row contains the lower errors, the second row contains the upper errors.
            # The reason xerr seems to be wrong way is that 'upper' refers to the upper ENERGY boundary, which corresponds to the LOWER 1/beta boundary
            if Onset and len(inverse_beta1)>0:
                label1 = Onset.sensor.upper() if mass_energy==mass_energy1 else f"{Onset.sensor.upper()} {species_title1}"
                ax.errorbar(inverse_beta1, onset_times1, yerr=y_errors1_plot, xerr=[x_errors_upper1, x_errors_lower1], 
                        fmt='o', elinewidth=VDA_ELINEWIDTH, capsize=VDA_CAPSIZE, zorder=1, label=label1)

            if not Onset:
                label = "onset times"
            else:
                label = self.sensor.upper() if mass_energy==mass_energy1 else f"{self.sensor.upper()} {species_title}"

            if len(inverse_beta) > 0:
                ax.errorbar(inverse_beta, onset_times, yerr=y_errors_plot, xerr=[x_errors_upper, x_errors_lower], 
                            fmt='o', elinewidth=VDA_ELINEWIDTH, capsize=VDA_CAPSIZE, zorder=1, label=label)

            # Omitted datapoints, paint all points white and then those not omitted blue (+ red) again
            if omitted_exists:
                if Onset:
                    ax.scatter(inverse_beta1[selection1], onset_times1[selection1], s=VDA_MARKERSIZE, zorder=3)

                ax.scatter(inverse_beta_all, onset_times_all, c="white", s=VDA_MARKERSIZE-5, zorder=2)
                ax.scatter(inverse_beta[selection], onset_times[selection], s=VDA_MARKERSIZE, zorder=3)

            # The odr fit
            # Here we need to first take the selection of i_beta_all and ONLY after that take the compressed form, which is the set of valid values
            ax.plot(inverse_beta_all[selection_all].compressed(), odr_fit, ls='--',
                label=f"L: {np.round(slope,3):.3f} $\pm$ {np.round(slope_uncertainty,3):.3f} AU\nt_inj: {rel_time_str} $\pm$ {str(t_inj_uncertainty)[7:]}")

            ax.set_xlabel(r"1/$\beta$", fontsize = AXLABEL_FONTSIZE)

            # Format the y-axis. For that make a selection to exclude NaTs from the set of onset times that define 
            # the vertical axis boundaries.
            nat_onsets = pd.isnull(onset_times_all)
            not_nats = np.array(onset_times_all)[~nat_onsets]

            # We have to catch on the above line all non-NaT onset times, because numpy nanmin() and nanmax() don't recognize them
            if len(not_nats) > 0:
                if np.nanmax(not_nats)-np.nanmin(not_nats) > pd.Timedelta(minutes=10):
                    y_axis_time_format = DateFormatter("%H:%M")
                    ax.set_ylabel("Time (HH:MM)", fontsize = AXLABEL_FONTSIZE)
                else:
                    y_axis_time_format = DateFormatter("%H:%M:%S")
                    ax.set_ylabel("Time (HH:MM:SS)", fontsize = AXLABEL_FONTSIZE)
            ax.yaxis.set_major_formatter(y_axis_time_format)

            if ylim:
                ax.set_ylim(pd.to_datetime(ylim[0]),pd.to_datetime(ylim[1]))

            set_standard_ticks(ax=ax)
            set_legend(ax=ax, legend_loc="in", fontsize=LEGEND_SIZE)

            # Title for the figure
            if title is None:
                if Onset:
                    instrument_species_id = instrument.upper() if species_title==species_title1 else f"{instrument.upper()} {species_title}"
                    # This is the default for joint VDA, two instruments of the same spacecraft
                    if self.spacecraft == Onset.spacecraft:
                        if self.viewing:
                            ax.set_title(f"VDA, {spacecraft.upper()} / {instrument_species_id} ({self.viewing}) + {Onset.sensor.upper()} {species_title1}, {date_of_event}", fontsize=TITLE_FONTSIZE)
                        else:
                            ax.set_title(f"VDA, {spacecraft.upper()} / {instrument_species_id} + {Onset.sensor.upper()} {species_title1}, {date_of_event}", fontsize=TITLE_FONTSIZE)

                    else:
                        # In this case these are two different spacecraft
                        if self.viewing and Onset.viewing:
                            ax.set_title(f"VDA, {spacecraft.upper()}/{instrument_species_id} ({self.viewing}) + {Onset.spacecraft.upper()}/{Onset.sensor.upper()}({Onset.viewing})\n{species_title1}, {date_of_event}", fontsize=TITLE_FONTSIZE)
                        elif self.viewing:
                            ax.set_title(f"VDA, {spacecraft.upper()}/{self.sensor} ({self.viewing}) {species_title} + {Onset.spacecraft.upper()}/{Onset.sensor.upper()}\n{species_title1}, {date_of_event}", fontsize=TITLE_FONTSIZE)
                        elif Onset.viewing:
                            ax.set_title(f"VDA, {spacecraft.upper()}/{instrument_species_id} + {Onset.spacecraft.upper()}/{Onset.sensor.upper()}({Onset.viewing})\n{species_title1}, {date_of_event}", fontsize=TITLE_FONTSIZE)
                        else:
                            ax.set_title(f"VDA, {spacecraft.upper()}/{instrument_species_id} + {Onset.spacecraft.upper()}/{Onset.sensor.upper()} {species_title1}, {date_of_event}", fontsize=TITLE_FONTSIZE)

                else:
                    # Single spacecraft, single instrument
                    if self.viewing:
                        ax.set_title(f"VDA, {spacecraft.upper()} / {instrument.upper()} ({self.viewing}) {species_title}, {date_of_event}", fontsize=TITLE_FONTSIZE)
                    else:
                        ax.set_title(f"VDA, {spacecraft.upper()} / {instrument.upper()} {species_title}, {date_of_event}", fontsize=TITLE_FONTSIZE)
            else:
                ax.set_title(title, fontsize=TITLE_FONTSIZE)

            # Saving the figure in "/plots"
            if save:
                if not savepath:
                    savepath = CURRENT_PATH
                if Onset:
                    savestr = f"{savepath}{os.sep}VDA_{spacecraft}_{instrument}({self.viewing})+{Onset.sensor}_{species_title}_{date_of_event}.png" if self.viewing else f"{savepath}{os.sep}VDA_{spacecraft}_{instrument}+{Onset.sensor}_{species_title}_{date_of_event}.png"
                else:
                    savestr = f"{savepath}{os.sep}VDA_{spacecraft}_{instrument}_{self.viewing}_{species_title}_{date_of_event}.png" if self.viewing else f"{savepath}{os.sep}VDA_{spacecraft}_{instrument}_{species_title}_{date_of_event}.png"
                plt.savefig(savestr, transparent=False,
                        facecolor="white", bbox_inches="tight")

        # Add the figure and axes of it to the return
        output["fig"] = fig
        output["axes"] = ax

        return output


    def automatic_onset_stats(self, channels, background, viewing, erase, cusum_minutes:int=None, sample_size:float=0.5, 
                              small_windows=None, stop=None, weights="inverse_variance", limit_computation_time=True, sigma=2, 
                              detrend:bool=True, prints:bool=False, custom_data_dt:str=None,
                              limit_averaging:str=None, fail_avg_stop:int=None):
        """
        Automates the uncertainty estimation for a single channel provided by Poisson-CUSUM-bootstrap hybrid method

        Parameters:
        ----------
        channels : {int}
                    Channel number.
        background : {BootstrapWindow}
                    The background that precedes the onset of the event.
        viewing : {str}
                    Viewing direction of the instrument.
        erase : {tuple(float, str)}
                    Acts a a filter that discards all values above the first entry BEFORE the timestamp defined by the second entry.
        cusum_minutes : {int}, optional
                    The amount of time in minutes that we demand on threshold-exceeding intensity.
        sample_size : {float} default 0.5
                    The fraction of the datapoints that the method will randomly pick.
        small_windows : {str} default None
                     Choose either 'random' or 'equal' to partition the big background window into smaller ones. 'equal' small windows
                     are equally placed inside the background, while 'random' placement allows the windows to overlap.
        stop : {str} default None
                    A pandas-compatible time string that decides on how long to continue averaging data while making onset 
                    distributions. If None, data will be averaged up to the 1-sigma uncertainty of whatever the native data resolution
                    produces.
        weights : {str} optional (default: 'inverse_variance')
                    Choose weights for calculating the mean uncertainty. 'inverse_variance' uses the inverse of the variance of the onset time distribution,
                    'int_time' will weight uncertainties by their individual integrating times, while 'uncertainty' will use the inverse of the ~95% of the
                    distribution as the weights.
        limit_computation_time : {bool}, default True
                    If enabled, skips all integration times above 10 minutes that are not multiples of 5. 
        sigma : {int, float} default 2
                    The multiplier for the $\mu_{d}$ variable in the CUSUM method.
        detrend : {bool}, default True
                    Switch to apply a shift on all but the native data distributions such that the onset times are shifted backwards in
                    time by half of the data time resolution.
        prints : {bool}, optional
                    Switch to print information about which channels is being analyzed and what's its 1-sigma uncertainty.
        custom_data_dt : {str}, optional
                    In case of custom data, provide the native cadence of the data PRIOR to resampling.
        limit_averaging : {str}, optional
                    Pandas-compatible time string. Limits the averaging to a certain time. Leave to None to not limit averaging. 
        fail_avg_stop : {int}, optional
                    If absolutely no onset is found in the native time resolution, how far should the method average the data to
                    try and find onset times? Default is up to 5 minutes.

        Returns:
        ----------
        stats_arr : {OnsetStatsArray}
        """

        self.background = background

        # A fine cadence is less than 1 minute -> requires computation time
        FINE_CADENCE_SC = ("solo", "wind")

        def dt_str_to_int(stop):
            """
            Converts stop condition string to a integer.
            """

            if isinstance(stop, str):
                    if stop[-3:] == "min" or stop[-1:] == 'T':
                        split_str = "min"
                        divisor = 1
                    elif stop[-1:] == 'H':
                        split_str = "h"
                        divisor = 1/60
                    elif stop[-1:] == 's':
                        split_str = 's'
                        divisor = 60
                    else:
                        raise ValueError(f"{stop} is an incorrect value for parameter 'average_to'")
                    stop_int = int(stop.split(split_str)[0])//divisor
            else:
                raise ValueError("Time string has to pandas-compatible time string, e.g., '15 min' or '60 s'.")
            
            return stop_int


        def produce_integration_times(int_time_ints, limit_averaging, stop) -> list[str]:
            """
            Returns a list of integration time strings to be used in the loop that produces
            onset statistics arrays.
            """

            if limit_averaging is not None:
                unit_of_time = find_biggest_nonzero_unit(timedelta=pd.Timedelta(limit_averaging))
            elif stop is not None:
                unit_of_time = find_biggest_nonzero_unit(timedelta=pd.Timedelta(stop))
            else:
                raise Exception("""The unit of integration time not identified! The issue could be caused \
                                by custom data with unidentified frequency. Try either limiting the averaging with 'limit_averaging' \
                                keyword, or averaging to a set frequency with 'average_to' keyword.""")
            
            int_time_strs = [f"{i} {unit_of_time}" for i in int_time_ints]

            return int_time_strs


        # Check channel validity (if standard data)
        if not self.custom_data:
            if not isinstance(channels,int):
                if not isinstance(channels,list):
                    channels = int(channels)

        # Get the integer numbers that represent stop and/or limit_averaging in minutes
        if stop:
            # Integer value of stopping condition, e.g., the '5' in '5 min'.
            stop_int = dt_str_to_int(stop)
            # Unit of time
            unit_of_time = find_biggest_nonzero_unit(timedelta=pd.Timedelta(stop))
        if limit_averaging:
            limit_averaging_int = dt_str_to_int(limit_averaging)
            # Unit of time
            unit_of_time = find_biggest_nonzero_unit(timedelta=pd.Timedelta(limit_averaging))

        # SolO/EPT first channel does not provide proper data as of late
        if self.spacecraft=="solo" and self.sensor=="ept" and channels==0:
            self.input_nat_onset_stats(channels)
            return  None

        # Wind/3DP first electron channel and the first two proton channels don't provide proper data
        if self.spacecraft=="wind" and self.species=='e' and channels==0:
            return  None
        if self.spacecraft=="wind" and self.species=='p' and channels in (0,1):
            return  None

        # The first round of onset statistics is acquired from 1 minute resolution, if computation time is limited 
        if self.spacecraft in FINE_CADENCE_SC and limit_computation_time:
            first_resample = "1 min"
        # Check if cadence is fine also for custom data
        elif self.custom_data:
            if custom_data_dt is None:
                custom_data_dt = self.data.index.freq if self.data.index.freq is not None else get_time_reso(self.data)

            # A fine cadence means less than 1 minute
            freq_is_fine = (pd.Timedelta(custom_data_dt) < pd.Timedelta("1 min"))

            first_resample = "1 min" if freq_is_fine and limit_computation_time else None

        else:
            first_resample = None

        # SOHO/EPHIN E300 is deactivated from 2017 onward -> there will be no reasonable onset there
        if self.spacecraft=="soho" and self.sensor=="ephin" and channels==300 and self.flux_series.index[0] >= EPHIN_300_INVALID_ONWARDS:
            print("Channel deactivated because of failure mode D.")
            return None

        # Run statistic_onset() once to get the confidence intervals for the bare not resampled, or 1-minute, data
        first_run_stats, _ = self.statistic_onset(channels=channels, Window=background, viewing=viewing, 
                                            sample_size=sample_size, resample=first_resample, erase=erase, small_windows=small_windows,
                                            cusum_minutes=cusum_minutes, detrend=False, sigma_multiplier=sigma)

        # For the first iteration initialize the OnsetStatsArray object, which can plot the integration time plot
        # This step has to be done after running statistic_onset() the first time, because otherwise "self.bootstrap_onset_statistics"
        # does not exist yet
        stats_arr = OnsetStatsArray(self)

        # Integer number of first run uncertainty in minutes
        first_run_uncertainty_mins = (first_run_stats["1-sigma_confidence_interval"][1] - first_run_stats["1-sigma_confidence_interval"][0]).seconds // 60
        first_run_uncertainty = first_run_stats["1-sigma_confidence_interval"][1] - first_run_stats["1-sigma_confidence_interval"][0]

        if prints:
            print(f"~68 % uncertainty for the self time with native data resolution: {first_run_uncertainty}")

        # Could be that the native resolution identifies no onset at all, in this case handle it
        if not isinstance(first_run_uncertainty, pd._libs.tslibs.nattype.NaTType):

            # If stop is not defined, then average up to predefined (default or first_run_uncertainty_mins) time
            if not stop:

                # Most of the high-energy particle instruments have a time resolution of 1 min, so don't do averaging for them
                # if uncertainty is something like 1 min 07 sec
                if first_run_uncertainty_mins < 2 and self.spacecraft not in FINE_CADENCE_SC:

                    stats_arr.calculate_weighted_uncertainty()
                    return stats_arr

                else:
                    # SolO instruments and Wind/3DP have high cadence (< 1 min), so start integrating from 1 minute measurements
                    # unless limit_computation_time is enabled
                    if not self.custom_data:
                        start_idx = 1 if (self.spacecraft in FINE_CADENCE_SC and not limit_computation_time) else 2

                    else:
                        # In case of custom data, start from the base cadence. Except for if the cadence is fine, then
                        # start from 1 minute
                        if freq_is_fine and not limit_computation_time:
                            start_idx = 1
                        
                        # If the cadence is not fine, then no matter if computation time is limited or not, we start
                        # time-averaging from the cadence_in_minutes+1 onward, e.g., 10 min -> 10+1 == 11 
                        else:
                            start_idx = pd.Timedelta(custom_data_dt).seconds//60 + 1

                    # Initialize integration times only up to the amount of minutes that the first run had uncertainty
                    int_times = np.array([i for i in range(start_idx,first_run_uncertainty_mins+1)], dtype=int)

                    if prints:
                        if limit_averaging:
                            upto_averaging_display = limit_averaging_int if limit_averaging_int < first_run_uncertainty_mins else first_run_uncertainty_mins
                        else:
                            upto_averaging_display = first_run_uncertainty_mins
                        print(f"Averaging up to {upto_averaging_display} minutes") if upto_averaging_display > 0 else print("Not averaging.")

            else:

                if stop_int > 0:
                    if prints:
                        print(f"Averaging up to {stop_int} minutes")
                    else:
                        pass

                # SolO instruments and Wind/3DP have high cadence (< 1 min), so start integrating from 1 minute measurements
                    # unless limit_computation_time is enabled
                if self.spacecraft in FINE_CADENCE_SC and not limit_computation_time:
                    int_times = np.array([i for i in range(1,stop_int+1)])
                else:
                    int_times = np.array([i for i in range(2,stop_int+1)])

        # Go here if no onset found at all
        else:

            if self.spacecraft in FINE_CADENCE_SC and not limit_computation_time:
                try_avg_start = 1
            else:
                try_avg_start = 2
            try_avg_stop = 5 if not isinstance(fail_avg_stop,int) else fail_avg_stop

            # Try up to {try_avg_stop} minutes averaging (default 5), if still no onset -> give up
            for i in range(try_avg_start,try_avg_stop+1):

                next_run_stats, _ = self.statistic_onset(channels=channels, Window=background, viewing=viewing, 
                                            sample_size=sample_size, resample=f"{i}min", erase=erase, small_windows=small_windows,
                                            cusum_minutes=cusum_minutes, sigma_multiplier=sigma, detrend=True)
                next_run_uncertainty = next_run_stats["1-sigma_confidence_interval"][1] - next_run_stats["1-sigma_confidence_interval"][0]
                next_run_uncertainty_mins = (next_run_stats["1-sigma_confidence_interval"][1] - next_run_stats["1-sigma_confidence_interval"][0]).seconds // 60

                if not isinstance(next_run_uncertainty, pd._libs.tslibs.nattype.NaTType):
                    if prints:
                        print(f"No onset found in the native data resolution. ~68 % uncertainty with {i} min resolution: {next_run_uncertainty}")

                    # Here check if it makes sense to average "from i minutes to <uncertainty> minutes or up to "stop" minutes
                    if stop:

                        if i < stop_int:
                            int_times = np.array([j for j in range(i,stop_int+1)])
                            if prints:
                                print(f"Averaging from {i} minutes up to {stop_int} minutes")
                        else:
                            if prints:
                                print(f"Stop condition set to {stop_int} minutes, which is less than {i} min. Using only {i} minutes averaged data.")
                            int_times = np.array([j for j in range(i,i+1)])

                    elif i < next_run_uncertainty_mins:
                        int_times = np.array([j for j in range(i,next_run_uncertainty_mins+1)])
                        if prints:
                            if limit_averaging:
                                upto_averaging_display = limit_averaging_int if limit_averaging_int < next_run_uncertainty_mins else next_run_uncertainty_mins
                            else:
                                upto_averaging_display = next_run_uncertainty_mins

                            print(f"Averaging from {i} minutes up to {upto_averaging_display} minutes")

                    else:
                        if prints:
                            print(f"~68 % uncertainty less than current time-averaging. Terminating.")

                        stats_arr.add(self)
                        stats_arr.calculate_weighted_uncertainty(weights)
                        return stats_arr

                    break

                else:
                    # If we tried everything and still no onset -> NaT and exit
                    if i==try_avg_stop:
                        if prints:
                            print(f"No onsets found with {i} min time averaging. Terminating.")
                            stats_arr.calculate_weighted_uncertainty("int_time")
                            return stats_arr
                    else:
                        pass


        # Here if int_times gets too long, coarsen it up a little from 15 minutes onward
        # Logic of this selection: preserve an element of int_times if it's at most 10 OR if it's divisible by 5
        if limit_computation_time:
            int_times = int_times[np.where((int_times <= 15) | (int_times%5==0))]

        # If the user set some upper limit to the averaging, apply that limit here
        if isinstance(limit_averaging,str):
            int_times = int_times[np.where(int_times <= limit_averaging_int)]

        # Finally convert int_times (integers) to pandas-compatible time strs
        int_time_strs = produce_integration_times(int_time_ints=int_times, limit_averaging=limit_averaging, stop=stop)           

        # Loop through int_times as far as the first run uncertainty reaches
        for resample in int_time_strs:

            # A check to adjust the cusum window in cases of large integration times. Only for standardized data.
            if not self.custom_data:
                # Here it is assumed that the resample string represents minutes, and ends with 'min'.
                if int(resample[:-3]) > 10:
                    cusum_minutes = int(resample[:-3])*3

            _, _ = self.statistic_onset(channels=channels, Window=background, viewing=viewing, 
                                            sample_size=sample_size, resample=resample, erase=erase, small_windows=small_windows,
                                            cusum_minutes=cusum_minutes, sigma_multiplier=sigma, detrend=detrend)

            stats_arr.add(self)

        # Calculate the weighted medians and confidence intervals. This method automatically updates the onset object's
        # dictionary of uncertainties as well.
        stats_arr.calculate_weighted_uncertainty(weights)

        return stats_arr


    def onset_statistics_per_channel(self, background, viewing, channels=None, erase:list=None, cusum_minutes:int=30, sample_size:float=0.50, 
                                     weights:str="inverse_variance", detrend=True, limit_computation_time=True, average_to=None, print_output=False, 
                                     limit_averaging=None, fail_avg_stop:int=None, random_seed:int=None, sigma:int=2):
        """
        Wrapper method for automatic_onset_stats(), that completes full onset and uncertainty analysis for a single channel.
        Does a complete onset uncertainty analysis on, by default all, the energy channels for the given instrument.

        Parameters:
        -----------
        channels : {str or array-like}, optional
                    A tuple, list or Range of channel ID's to run the method on. Leave to None or set to 'all' to run for all channels.
        background : {BootstrapWindow}
                    The common pre-event background used for the energy channels. 
        viewing : {str}
                    The viewing direction if the instrument.
        erase : {tuple(float, str)}, optional
                    If there are spikes in the background one wishes to omit, set the threshold value and ending point for ignoring those
                    values here. Example: [1000, '2022-05-10 12:00'] discards all values 1000 or above until 2022-05-10 noon.
        cusum_minutes : {int, float}, optional
                    The amount of MINUTES that the method will demand for continuous threshold-axceeding intensity before identifying
                    an onset. 
        sample_size : {float}, optional
                    The fraction of the data points inside the background window that will be considered for each of the bootstrapped
                    runs of the method.
        sigma : {int}, optional
                    The multiplier for the standard deviation, sigma, in the Poisson-CUSUM method. Defaults to 2.
        weights : {str}, optional (default: 'inverse_variance')
                    Pick 'inverse_variance' to use inverse variance weighting, 'uncertainty' to use the width of 2-sigma intervals as 
                    the basis for weighting timestamps, or 'int_time' to use the integration time as a basis for the weighting.
        detrend : {bool}, optional
                    Switch to apply detrending on the onset time distributions.
        limit_computation_time : {bool}, optional
                    If True (default), then resample all time series that have finer time resolution than 1 minute, to 1 minute resolution
                    before 
        average_to : {str}, optional
                    Explicitly tells the method to average every channel up to a specific time resolution, disregarding the
                    recommendation got from 1-sigma width of the native data. If both 'average_to' and 'limit_averaging' are
                    given as an input, 'limit_averaging' will take precedence over 'average_to'.
        print_output : {bool}
                    Switch to print when a new channel is being analyzed and for how far it will be time-averaged to.
        limit_averaging : {str}, optional
                    Pandas compatible time string to limit the averaging time to a certain time, e.g., '60 min'
        fail_avg_stop : {int}, optional
                    If absolutely no onset is found in the native time resolution, how far should the method average the data to
                    try and find onset times? Default is up to 5 minutes.
        random_seed : {int}, optional
                    Passes down a seed for the random generator that picks the samples from the background window.

        Returns:
        ----------
        uncertainty_stats_by_channel : {np.ndarray(OnsetStatsArray)}
                    A numpy array of OnsetStatsArray objects, each of which encapsulates statistics of the onset wime in each of the channels
                    that the method was run over. 
        """

        # If a random seed (and a valid one!) was given, apply it before doing anything else
        if isinstance(random_seed,int):
            np.random.seed(random_seed)

        # Initialize the array which will store all the uncertaintystats objects
        uncertainty_stats_by_channel = np.array([])

        # Check which channels to run uncertainty stats on
        if channels is None or channels=="all":
            if self.custom_data:
                raise Exception("The list of channels must be manually given as an input in the case of custom data!")
            all_channels = self.get_all_channels()
        elif isinstance(channels, (tuple, list, range)):
            all_channels = channels
        elif isinstance(channels, (int, np.int64)):
            all_channels = [channels]
        elif isinstance(channels,str) and self.custom_data:
            all_channels = [channels]
        else:
            raise TypeError(f"{type(channels)} is and incorrect type of argument 'channels'! It should be None, str=='all', tuple, list or range.")

        # Recognize the base cadence here, before going in to automatic_onset_stats
        if self.custom_data:

            # Either the user hasn't provided a frequency to the data, or it's not possible due to irregular 
            # timedelta between data points. Try to infer the most prevalent timedelta, and assume that it is
            # a good representation of the time resolution.
            if self.data.index.freq is None:
                custom_data_dt = get_time_reso(self.data)

                if print_output:
                    print("Data frequency doesn't appear to be set. This can be done with set_data_frequency() -method.")
                    print(f"Inferring the most probable time resolution: {custom_data_dt}")
            else:
                custom_data_dt = self.data.index.freq
        else:
            custom_data_dt = None

        if print_output:
            background.print_max_recommended_reso()

        # Loop through all the channels and save the onset statistics to objects that will be stored in the array initialized in the start
        for channel in all_channels:

            if print_output:
                print(f"Channel {channel}:")

            # automatic_onset_stats() -method runs statistic_onset() for all different data integration times for a single channel
            onset_uncertainty_stats = self.automatic_onset_stats(channels=channel, background=background, viewing=viewing, erase=erase, sigma=sigma,
                                                                stop=average_to, cusum_minutes=cusum_minutes, sample_size=sample_size, 
                                                                weights=weights, detrend=detrend, limit_computation_time=limit_computation_time,
                                                                prints=print_output, limit_averaging=limit_averaging, fail_avg_stop=fail_avg_stop,
                                                                custom_data_dt=custom_data_dt)

            # Add statistics to the array that holds statistics related to each individual channel
            uncertainty_stats_by_channel = np.append(uncertainty_stats_by_channel, onset_uncertainty_stats)

        return uncertainty_stats_by_channel


    def find_first_peak(self, channels, window, xlim, viewing, resample, plot=True):
        """
        Finds the first maximum in time series data.
        
        Parameters:
        -----------
        window : {tuple(datetime str)}
        
        Returns:
        ----------
        first_maximum_peak_time : {pd.datetime}
        """

        flux_series, en_channel_string = self.choose_flux_series(channels=channels, viewing=viewing)

        # Resample data if requested
        if resample is not None:
            flux_series = util.resample_df(flux_series, resample)

        # Cut flux_series according to given window
        window = pd.to_datetime(window)
        cut_series = flux_series.loc[(flux_series.index > window[0]) & (flux_series.index < window[1])]

        # The series is indexed by time
        time = cut_series.index

        first_maximum_peak_val = cut_series.values[0]
        first_maximum_peak_time = time[0]
        for i, val in enumerate(cut_series.values):

            if val > first_maximum_peak_val:
                first_maximum_peak_val = val
                first_maximum_peak_time = time[i]

        # Plotting commands ->
        fig, ax = plt.subplots(figsize=(21,9))

        # Setting the x-axis limits
        if xlim is None:
            xlim = [flux_series.index[0], flux_series.index[-1]]
        else:
            # Check that xlim makes sense
            if xlim[0] == xlim[1]:
                raise Exception("xlim[0] and xlim[1] cannot be the same time!")

            xlim[0], xlim[1] = pd.to_datetime(xlim[0]), pd.to_datetime(xlim[1])
            ax.set_xlim(xlim)

        ax.set_yscale("log")

        # The measurements. kw where dictates where the step happens -> "mid" for at the middle of the bin
        ax.step(flux_series.index, flux_series.values, c="C0", where="mid")

        # Window shaded area
        ax.axvspan(xmin=window[0], xmax=window[1],
                    color="lightblue", label="Window", alpha=0.3)

        # Peak time, value and a textbox displaying both
        ax.axvline(x=first_maximum_peak_time, linewidth=1.2, color="navy", linestyle= '-', 
                    label="First maximum peak")

        # Textbox settings:
        plabel = AnchoredText(f"First maximum peak time+value\n{first_maximum_peak_time}\n{first_maximum_peak_val}", prop=dict(size=13), frameon=True, loc=(4) )
        plabel.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        plabel.patch.set_linewidth(2.0)
        ax.add_artist(plabel)

        # Labels for axes
        ax.set_xlabel("Time", fontsize=AXLABEL_FONTSIZE)
        ax.set_ylabel(f"Intensity [1/(cm^2 sr s MeV)]", fontsize=AXLABEL_FONTSIZE)

        # Tickmarks, their size etc...
        ax.tick_params(which='major', length=6, width=2, labelsize=17)
        ax.tick_params(which='minor', length=5, width=1.4)

        # Date tick locator and formatter
        ax.xaxis_date()
        utc_dt_format1 = DateFormatter('%H:%M \n%Y-%m-%d')
        ax.xaxis.set_major_formatter(utc_dt_format1)

        # Setting the title
        if self.species in ["electron", 'e']:
            s_identifier = 'electrons'
        if self.species in ["proton", 'p', 'H']:
            s_identifier = 'protons'
        if self.species in ["ion", 'i']:
            s_identifier = 'ions'

        if (viewing != '' and viewing is not None):

            plt.title(f"{self.spacecraft}/{self.sensor.upper()} {en_channel_string} {s_identifier}\n"
                    f"{resample} averaging, viewing: "
                    f"{viewing.upper()}")

        else:

            plt.title(f"{self.spacecraft}/{self.sensor.upper()} {en_channel_string} {s_identifier}\n"
                    f"{resample} averaging")

        if plot:
            plt.show()
        else:
            plt.close()

        return first_maximum_peak_time


    def get_distribution_statistics(self, onset_statistics:dict, percentiles:tuple = [(15.89,84.10), (2.30, 97.7)]):
        """
        Returns the median, mode, mean and confidence intervals of an onset distribution.
        
        Parameters:
        -----------
        onset_statistics : dict
                        A bootstrap_onset_statistics dictionary produced by the Onset class
        percentiles : tuple or a list of tuples
                        A tuple of percentiles, or a list containing the pairs of percentiles.
        
        Returns:
        -----------
        median : pandas datetime
                The median of the onset distribution
        mode : pandas datetime
                The mode of the onset distribution
        mean : pandas datetime
                The mean of the onset distribution
        confidence_intervals : tuple or a list of tuples
                Contains the percentiles or a list of percentiles if input was a list
        """

        onsets = pd.to_datetime(onset_statistics["onset_list"])

        # Calculate median with numpy.median
        median =  onsets[len(onsets)//2]

        # The most frequently appearing onset is already found
        mode = onset_statistics["unique_onsets"][-1,0]

        # Mean onset is already calculated
        mean = onset_statistics["mean_onset"]

        # Find the percentiles that were given
        confidence_intervals = self.get_distribution_percentiles(onsets, percentiles)

        return median, mode, mean, confidence_intervals


    def get_distribution_percentiles(self, onset_list, percentiles:tuple, time_reso=None):
        """
        Returns the confidence intervals of an onset distribution. Includes a check to make sure
        that the uncertainty intervals are not smaller than the data time resolution.
        
        Parameters:
        -----------
        percentiles : tuple or a list of tuples
                        A tuple of percentiles, or a list containing the pairs of percentiles.
        
        Returns:
        -----------
        confidence_intervals : tuple or a list of tuples
                Contains the percentiles or a list of percentiles if input was a list
        """

        # Find the percentiles that were given
        if isinstance(percentiles[0],(tuple,list)):
            confidence_intervals = []
            for pair in percentiles:
                conf_interval = pd.to_datetime(np.nanpercentile(onset_list,pair))
                confidence_intervals.append(conf_interval)
        else:
            confidence_intervals = pd.to_datetime(np.nanpercentile(onset_list,percentiles))

        # Finally check that the percentiles make sense (in that they are not less than the native data resolution of the instrument)
        if not time_reso:
            time_reso = self.get_minimum_cadence()
        confidence_intervals = check_confidence_intervals(confidence_intervals, time_reso=time_reso)

        return confidence_intervals


    def tsa_per_channel(self, radial_distance=None, path_length=None, solar_wind_speed=400, onset_times=None) -> dict:
        """
        Applies a time shift to all available energy channels accoridng to their kinetic energies
        and an assumed path of Parker spiral arc that they travelled.

        Parameters:
        -----------
        radial_distance : {float} The radial distance of the observer from the Sun at the time of observation.

        path_length : {float} Directly give the path length -> no calculation needed. This variable takes 
                              precedence over radial_distance.

        solar_wind_speed : {int/float} Solar wind speed in km/s. If not given, default to 400 km/s.

        onset_times : {dict} Give input onset times as a dictionary, coupling the channel identifier to the timestamp. 
                             If not given, use self.onset_statistics. (default=None)

        Returns: 
        ----------
        tsa_times : {dict} The time-shifted timestamps.
        """

        if radial_distance is None and path_length is None:
            raise TypeError("Either 'radial_distance' or 'path_length' must be specified to apply TSA!")

        # The nominal Parker spiral length for a given radial distance and solar wind speed
        if path_length is None:
            path_length = path_length_calculator(distance=radial_distance, solar_wind=solar_wind_speed)

        # The mean speeds of the energetic particles. The problem here is that the array
        # requires knowledge of the channel indices.
        try:
            sep_speeds = self.calculate_particle_speeds()

        # The UnboundLocalError arises due to seppy.calculate_particle_speeds() not recognizing the
        # energy values of the channels. Most likely because user has some custom data.
        except UnboundLocalError:
            sep_speeds = self.calculate_particle_speeds_custom()

        if onset_times is None:
            onset_times = self.onset_statistics
        else:
            # Check that every timestamp is in a list
            try:
                _ = onset_times[list(onset_times.keys())[0]][0]
            except TypeError:
                onset_times = {channel:[timestamp] for (channel,timestamp) in onset_times.items()}

        # Init a dictionary to collect time-shifted timestamps
        tsa_times = {}

        # Enumerating a dictionary yields i, the index of a key, and the key itself.
        for i, channel in enumerate(onset_times):
            onset_stats = onset_times[channel]
            tsa_times[channel] = tsa(t0 = onset_stats[0], L=path_length, v=sep_speeds[i])

        return tsa_times


    def tsa_plot(self, radial_distance=None, path_length=None, solar_wind_speed=400, ylim=None, 
                 plot=True, save=False, savepath=None, onset_times=None):
        """
        Applies tsa on all channels, and produces a scatter plot.
        
        Parameters:
        -----------
        radial_distance : {float} Radial heliocentric distance in AUs. Either this or <path_length> has to be specified.

        path_length : {float} Distance for TSA in AUs. Either this or <radial_distance> has to be specified.
        
        solar_wind_speed : {int/float} Speed of the solar wind for Parker spiral caclulation in km/s.
                                        Defaults to 400 km/s if not given.
        
        ylim : {tuple/list} Set custom limits for y-axis.
        
        plot : {bool} Switch to produce a plot.
        
        save : {bool} Switch to save the plot.

        savepath : {str} A path to save the plot, if the parameter <save> is enabled.

        onset_times : {dict} User-input onset times. If not give, default to using self.onset_statistics
        
        Returns:
        --------
        tsa_results : {dict} 
        """

        if radial_distance is None and path_length is None:
            raise TypeError("Either 'radial_distance' or 'path_length' must be specified for TSA!")

        species_str = "electron" if self.species in ELECTRON_IDENTIFIERS else "proton"

        tsa_results = {}

        # Gets the time-shifted timestamps and save them to the dictionary
        tsa_timestamps = self.tsa_per_channel(radial_distance=radial_distance, path_length=path_length, 
                                              solar_wind_speed=solar_wind_speed, onset_times=onset_times)
        tsa_results["tsa_timestamps"] = tsa_timestamps

        # The x-axis in terms of the inverse speed:
        try:
            inverse_betas = np.array([const.c.value/v for v in self.calculate_particle_speeds()])

        # The UnboundLocalError arises due to seppy.calculate_particle_speeds() not recognizing the
        # energy values of the channels. Most likely because user has some custom data.
        except UnboundLocalError:
            inverse_betas = np.array([const.c.value/v for v in self.calculate_particle_speeds_custom()])

        # Stupid check and not general in its nature, but solo first channel is unavailable
        # so leave it out here
        if self.spacecraft.lower()=="solo" and self.sensor=="ept":
            inverse_betas = inverse_betas[1:]

        # Init the figure
        tsa_fig, tsa_ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        # tsa_timestamps is a dictionary, so remember to only get the values for plotting
        tsa_ax.scatter(inverse_betas, tsa_timestamps.values(), s=135)

        tsa_ax.set_title(f"{self.spacecraft.upper()} / {self.sensor.upper()} {species_str} TSA", fontsize=TITLE_FONTSIZE)

        tsa_ax.yaxis.set_major_formatter(DateFormatter("%H:%M"))

        set_standard_ticks(ax=tsa_ax)
        
        tsa_ax.set_ylabel("Time", fontsize=AXLABEL_FONTSIZE)
        tsa_ax.set_xlabel(r"1/$\beta$", fontsize=AXLABEL_FONTSIZE)
        
        if ylim:
            tsa_ax.set_ylim(pd.to_datetime(ylim[0]), pd.to_datetime(ylim[1]))
        
        if save:
            if not savepath:
                savepath = CURRENT_PATH
            tsa_fig.savefig(f"{self.spacecraft.lower()}_{self.sensor.lower()}_{self.species}_tsa.png", facecolor="white", transparent=False,
                        bbox_inches="tight")
        if plot:
            plt.show()
        else:
            plt.close()

        tsa_results["inverse_betas"] = inverse_betas
        tsa_results["fig"] = tsa_fig
        tsa_results["axes"] = tsa_ax

        return tsa_results

    def calculate_particle_speeds_custom(self):
        """
        Calculates the average particle speeds by user-input channel energy boundaries.

        Note that the method set_custom_channel_energies() must have been ran before this method can be utilized.
        """

        if self.species in ELECTRON_IDENTIFIERS:
            m_species = const.m_e.value
        elif self.species in PROTON_IDENTIFIERS:
            m_species = const.m_p.value
        else:
            raise ValueError(f"The particle species {self.species} does not appear to be any of the recognized particle species. Can not calculate particle energy.")

        # E = mc^2, a fundamental property of any object with mass
        mass_energy = m_species * C_SQUARED
        
        # Get the energies of each energy channel, to calculate the mean energy of particles and ultimately
        # To get the dimensionless speeds of the particles (beta)
        e_lows, e_highs = self.get_custom_channel_energies()

        mean_energies = np.sqrt(np.multiply(e_lows, e_highs))

        # Transform kinetic energy from electron volts to joules
        e_Joule = [((En*u.eV).to(u.J)).value for En in mean_energies]

        # Beta, the unitless speed (v/c)
        beta = [np.sqrt(1-((e_J/mass_energy + 1)**(-2))) for e_J in e_Joule]

        return np.array(beta)*const.c.value


# The class that holds background window length+position and bootstrapping parameters
class BootstrapWindow:

    def __init__(self, start, end, n_shifts=0, bootstraps=1):
        self.start = pd.to_datetime(start)
        self.end = pd.to_datetime(end)
        self.n_shifts = n_shifts
        self.bootstraps = bootstraps

        self.attrs_dict = {
                        "start" : self.start,
                        "end" : self.end,
                        "n_shifts" : self.n_shifts,
                        "bootstraps": self.bootstraps
                          }

        self.init_state = self.save_state(return_state=True)

    def __repr__(self):
        return self.attrs()

    def __len__(self):
        return int((self.end - self.start).total_seconds()//60)

    def attrs(self, key=None):

        if key is None:
            return str(self.attrs_dict)
        else:
            return self.attrs_dict[key]

    def draw_background(self, ax):
        """
        Draws the window on a plot, given the axes of a plot.

        Parameters:
        -----------
        ax : {plt.Axes}
                The axes of the figure.
        """
        ax.axvspan(xmin=self.start, xmax=self.end,
                        color="#e41a1c", label="Background", alpha=BACKGROUND_ALPHA)
    
        # Textbox indicating onset time
        blabel = AnchoredText(f"Background:\n{self.start} - {self.end}", prop=dict(size=TXTBOX_SIZE), frameon=True, loc="upper left")
        blabel.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        blabel.patch.set_linewidth(2.0)
        ax.add_artist(blabel)

    def print_max_recommended_reso(self):
        f"""
        Prints out the maximum recommended resolution that the time series should be averaged to in order to still have
        at least {MIN_RECOMMENDED_POINTS} data points inside the background window.
        """

        minutes_in_background = len(self)

        # We recommend a maximum reso such that there are at least {MIN_RECOMMENDED_POINTS} data points to pick from
        max_reso = int(minutes_in_background/MIN_RECOMMENDED_POINTS)
        print(f"Your chosen background is {minutes_in_background} minutes long. To preserve the minimum of {MIN_RECOMMENDED_POINTS} data points to choose from,\nit is recommended that you either limit averaging up to {max_reso} minutes or enlarge the background window.")


    def move_window(self, hours=0, minutes=0):
        """
        Moves the BootstrapWindow by a given amount of hours and minutes.
        """
        if self.n_shifts > 0:
            self.start = self.start + pd.Timedelta(hours=hours, minutes=minutes)
            self.end = self.end +pd.Timedelta(hours=hours, minutes=minutes)
            self.n_shifts = self.n_shifts - 1

    def save_state(self, return_state=False):
        """
        Saves the state in which the window object is at the time of running this method.
        Mainly used to save / load the initialization state

        Parameters:
        -----------
        return_state : {bool}, default False
                    if True, return the dictionary that contains the state of class attributes. This is used to save
                    the initialization state of the object to a class attribute.
        """

        self.saved_state = {
            "start" : self.start,
            "end" : self.end,
            "n_shifts" : self.n_shifts,
            "bootstraps" : self.bootstraps
            }

        if return_state:
            return self.saved_state

    def load_state(self, which="previous"):
        """
        Loads the by default the previously saved state, i.e., sets the object's attributes back to what they were at the previous save point.
        """

        if which == "previous":

            self.start = self.saved_state["start"]
            self.end = self.saved_state["end"]
            self.n_shifts = self.saved_state["n_shifts"]
            self.bootstraps = self.saved_state["bootstraps"]

        elif which == "init":

            self.start = self.init_state["start"]
            self.end = self.init_state["end"]
            self.n_shifts = self.init_state["n_shifts"]
            self.bootstraps = self.init_state["bootstraps"]

        else:

            raise Exception(f"Can only load either 'previous' or 'init' state! The input was {which}")



class OnsetStatsArray:
    """
    Contains statistics and uncertainty of a single onset at a particular energy channel.
    """

    def __init__(self, onset_object):

        # Initialize the list of statistics (NOT A LIST OF ONSET OBJECTS)
        # self.statistics always points to the most recent dictionary in the archive of dictionaries
        self.statistics = onset_object.bootstrap_onset_statistics.copy()
        self.archive = [self.statistics]

        # Save the individual time series for later plotting purposes
        self.list_of_series = [onset_object.flux_series]

        # Save the onset object attributes
        self.spacecraft = onset_object.spacecraft
        self.sensor = onset_object.sensor
        self.species = onset_object.species
        self.channel_str = onset_object.recently_examined_channel_str
        self.channel_id = onset_object.current_channel_id

        # Remember to which onset object this statistics array is linked to
        self.linked_object = onset_object

        # Init integration times as an empty list and immediately calculate the first integration time
        self.integration_times = []
        self.calculate_integration_time(onset_object=onset_object)


    def __repr__(self):
        return str(f"Dictionaries in archive: {len(self.archive)}")

    def __len__(self):
        return len(self.archive)

    def attrs(self):
        print("OnsetStatsArray.statistics: the most recently added dictionary in the archive of dictionaries")
        print("OnsetStatsArray.archive: a list containing all the added dictionaries that each contain onset statistics")
        print("OnsetStatsArray.integration_times: a list containing integration times <float> in minutes corresponding to the statistics added to the archive")


    def add(self, onset_object):
        """
        Adds the bootstrap statistics of an onset object into the 'archive' of OnsetStatsArray. In practicality the archive is just a 
        list of dictionaries.
        """

        # Check that these are all the same sc / instrument / particle species
        if len(self.archive) >= 1:
            if onset_object.spacecraft.lower() != self.spacecraft.lower():
                raise Exception("Only one spacecraft per OnsetStatsArray permitted!")
            if onset_object.sensor.lower() != self.sensor.lower():
                raise Exception("Only one sensor per OnsetStatsArray permitted!")
            if onset_object.species.lower() != self.species.lower():
                raise Exception("Only one particle species per OnsetStatsArray permitted!")

        # Assert statistics as the latest stats dictionary
        self.statistics = onset_object.bootstrap_onset_statistics.copy()

        # Add the latest statistics dictionary to the list of dictionaries
        self.archive.append(self.statistics)

        # Calculate the most recent integration time and add it to the list of integration times
        self.calculate_integration_time(onset_object=onset_object)

        # Add the latest flux_series to the list of series for plotting purposes
        self.list_of_series.append(onset_object.flux_series)


    def drop(self, index:int=-1):
        """
        Drops by default the latest dictionary from the archive, unless given an index, 
        in which case the onset corresponding to the index is dropped. Also gets rid of 
        the chosen flux_series from the list of series.
        """
        _ = self.archive.pop(index)
        _ = self.integration_times.pop(index)
        _ = self.list_of_series.pop(index)

    def set_mean_of_median_onsets(self, mean_of_medians_onset):

        self.mean_of_medians_onset = mean_of_medians_onset

        self.linked_object.mean_of_medians_onset_acquired = True
        self.linked_object.mean_of_medians_onset = mean_of_medians_onset

    def set_w_median_and_confidence_intervals(self, mode, median, conf1_low, conf1_high, conf2_low, conf2_high):
        """
        Sets the weighted median and confidence intervals for onset time 
        """

        # Save the median and class attributes to the linked object's database:
        self.linked_object.update_onset_statistics([mode, median, conf1_low, conf1_high, conf2_low, conf2_high])

    def calculate_integration_time(self, onset_object) -> None:
        """
        Calculates the integration time (i.e. time resolution of data) of the given onset object and adds it to the list of integration times.
        Saves the integration times to a class attribute integration_times, which is a list of integration times in minutes as floats.
        """
        delta_sec = (onset_object.flux_series.index[2] - onset_object.flux_series.index[1]).seconds
        delta_min = delta_sec/60
        int_time = f"{delta_min} min" if delta_sec > 59 else f"{delta_sec} s"
        self.integration_times.append(int_time) if str(delta_min)[:7] != "0.98333" else self.integration_times.append("1 min")


    def onset_time_histogram(self, integration_time_index=0, binwidth="1 min", xlim=None, ylims=None, legend_loc=1,
                            save=False, savepath=None, grid=True):
        """
        A method to display the probability density histogram for the distribution of onset times collected
        to the object.

        Parameters:
        -----------
        integration_time_index : {int}, default 0
                                    Chooses the integration time, where 0 represents the native data time resolution.
        binwidth : {str}, default '1 min'
                                    Sets the width of time bins on the x-axis.
        xlims : {tuple}
        ylims : {tuple}
        legend_loc : {int}, default 1
                                    Sets the location of the legend.
        save : {bool}, default False
                                    Boolean save switch.
        savepath : {str}, optional
                                    Path to the directory to save the figure.
        grid : {bool}, default True
                                    Boolean switch to apply gridlines.
        """

        stats = self.archive[integration_time_index]

        # Plotting 
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        # If xlims not manually defined, let them be \pm 2 minutes from the first and last onset of the distribution
        if not xlim:
            # There might be NaTs among the onset list, so here choose only the valid ones to avert 
            # TypeError when taking np.nanmin() of the list
            valid_onset_list = np.array(stats["onset_list"])[~pd.isnull(stats["onset_list"])]
            xlim = (np.nanmin(valid_onset_list) - pd.Timedelta(minutes=2), np.nanmax(valid_onset_list) + pd.Timedelta(minutes=2))
        ax.set_xlim(xlim)

        # Show percentage on y-axis
        yvalues = [m/10 for m in range(0,11)]
        ax.set_yticks(yvalues)
        ax.set_yticklabels(['{}'.format(np.round(x, 1)) for x in yvalues])


        # Bins for the x axis 
        half_bin = pd.Timedelta(seconds=30)
        bins = pd.date_range(start=xlim[0]+half_bin, end=xlim[1]+half_bin, freq=binwidth).tolist()

        onset_frequencies = np.ones_like(stats["onset_list"])/len(stats["onset_list"])

        # Plotting the histogram rectangles
        bar_heights, bins, patches = ax.hist(stats["onset_list"], bins=bins, color="lightblue", edgecolor="black", weights=onset_frequencies, zorder=2)

        if not ylims:
            ylims = (0, np.nanmax(bar_heights)+0.02)
        ax.set_ylim(ylims)
    
        set_standard_ticks(ax=ax)

        # Mean, mode and median onset times as vertical solid lines
        ax.axvline(stats["mean_onset"], linewidth=2.0, color=COLOR_SCHEME["mean"], label=f"mean {str(stats['mean_onset'].time())[:8]}", zorder=3)
        ax.axvline(stats["most_likely_onset"][0], linewidth=2.0, color=COLOR_SCHEME["mode"], label=f"mode {str(stats['most_likely_onset'][0].time())[:8]}", zorder=3)
        ax.axvline(stats["median_onset"], linewidth=2.0, color=COLOR_SCHEME["median"], label=f"median {str(stats['median_onset'].time())[:8]}", zorder=3)

        # 1 -and 2-sigma intervals as red and blue dashed lines
        ax.axvspan(xmin=stats["2-sigma_confidence_interval"][0], xmax=stats["2-sigma_confidence_interval"][1], color=COLOR_SCHEME["2-sigma"], alpha=0.15, label="~95 % confidence", zorder=1)
        ax.axvspan(xmin=stats["1-sigma_confidence_interval"][0], xmax=stats["1-sigma_confidence_interval"][1], color=COLOR_SCHEME["1-sigma"], alpha=0.15, label="~68 % confidence", zorder=1)

        ax.set_xlabel(f"Time ({stats['mean_onset'].strftime('%Y-%m-%d')})", fontsize=AXLABEL_FONTSIZE)
        ax.set_ylabel("PD", fontsize=AXLABEL_FONTSIZE)

        ax.grid(visible=grid, axis="both")

        integration_time_str = f"{self.integration_times[integration_time_index]} integration time" if pd.Timedelta(self.integration_times[integration_time_index]) != self.linked_object.get_minimum_cadence() else f"{self.integration_times[integration_time_index]} data"

        ax.set_title(f"Probability density for {self.linked_object.background.bootstraps} onset times\n{integration_time_str}", fontsize=TITLE_FONTSIZE)
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))

        ax.legend(loc=legend_loc, bbox_to_anchor=(1.0, 1.0), fancybox=True, ncol=3, fontsize = LEGEND_SIZE)
        #set_legend(ax=ax, legend_loc="in", fontsize=LEGEND_SIZE)

        if save:
            if savepath is None:
                plt.savefig("onset_times_histogram.png", transparent=False, facecolor="white", bbox_inches="tight")
            else:
                plt.savefig(f"{savepath}{os.sep}onset_times_histogram.png", transparent=False, facecolor="white", bbox_inches="tight")

        plt.show()

    def integration_time_plot(self, title=None, ylims=None, save=False, savepath:str=None, grid=True) -> None:
        """
        Plots the median, mean, mode and confidence intervals for an array of Onset objects as a function
        of integration time (basically time resolution of the data, may be resampled)

        Parameters:
        -----------
        title : {str}, default None
                Title of the figure, generate automatically if None
        ylims : {array-like}, len==2, optional
                The lower and upper limit of the y-axis as datetimes or pandas-compatible datetime strings,
                e.g., '2010-12-31 14:00'.
        save : {bool}, default False
                Saves the figure
        savepath : {str}, optional
                The directory path or subdirectory to save the figure.
        grid : {bool}, default True
                Boolean switch to turn on gridlines.
        """

        # Collect the stats and different time resolutions
        means = []
        medians = []
        modes = []
        confidence_interval_1sigma = []
        confidence_interval_2sigma = []
        for statistics in self.archive:

            means.append(statistics["mean_onset"])
            medians.append(statistics["median_onset"])
            modes.append(statistics["most_likely_onset"][0])
            confidence_interval_1sigma.append(statistics["1-sigma_confidence_interval"])
            confidence_interval_2sigma.append(statistics["2-sigma_confidence_interval"])

        # This is done to preserve the confidence_intervals as class attributes that can be used later
        self.confidence_intervals_1sigma = np.array(confidence_interval_1sigma)
        self.confidence_intervals_2sigma = np.array(confidence_interval_2sigma)

        # Lower and higher boundaries of the confidence intervals of the distributions, in order of appearance
        conf1_lows = [pair[0] for pair in confidence_interval_1sigma]
        conf1_highs = [pair[1] for pair in confidence_interval_1sigma]

        conf2_lows = [pair[0] for pair in confidence_interval_2sigma]
        conf2_highs = [pair[1] for pair in confidence_interval_2sigma]

        # these are the real weighted mode and median
        mean_of_modes = self.linked_object.onset_statistics[self.channel_id][0]
        mean_of_medians = self.linked_object.onset_statistics[self.channel_id][1]

        figdate = get_figdate(modes)

        # Initializing the figure
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        # These integration time strings are transformed to floating point numbers in units of minutes
        xaxis_int_times = [pd.Timedelta(td).seconds/60 for td in self.integration_times]

        ax.set_xlim(0, xaxis_int_times[-1]+1)

        if ylims:
            ax.set_ylim(pd.to_datetime(ylims[0]), pd.to_datetime(ylims[1]))

        # We only want to have integer tickmarks
        ax.set_xticks(range(0,int(xaxis_int_times[-1]+1)))

        ax.scatter(xaxis_int_times, means, s=165, label="mean", zorder=2, color=COLOR_SCHEME["mean"], marker=".")
        ax.scatter(xaxis_int_times, medians, s=115, label="median", zorder=2, color=COLOR_SCHEME["median"], marker="^")
        ax.scatter(xaxis_int_times, modes, s=115, label="mode", zorder=2, color=COLOR_SCHEME["mode"], marker="p")

        ax.axhline(y=mean_of_medians, color=COLOR_SCHEME["median"], lw=2, label=f"Mean of medians:\n{str(mean_of_medians.time())[:8]}")
        ax.axhline(y=mean_of_modes, color=COLOR_SCHEME["mode"], lw=2, label=f"Mean of modes:\n{str(mean_of_modes.time())[:8]}")

        ax.fill_between(xaxis_int_times, y1=conf1_lows, y2=conf1_highs, facecolor=COLOR_SCHEME["1-sigma"], alpha=0.3, zorder=1)
        ax.fill_between(xaxis_int_times, y1=conf2_lows, y2=conf2_highs, facecolor=COLOR_SCHEME["2-sigma"], alpha=0.3, zorder=1)

        # For some reason all the labelsizes appear smaller in the plots created with this method, which is 
        # why I add some extra size to them.

        ax.set_xlabel("Data integration time [min]", fontsize=AXLABEL_FONTSIZE+6)
        ax.set_ylabel(f"{figdate}\nTime [HH:MM]", fontsize=AXLABEL_FONTSIZE+6)

        if not title:
            # particle_str = "electrons" if self.species=='e' else "protons" if self.species=='p' else "ions"
            ax.set_title(f"{self.spacecraft.upper()} / {self.sensor.upper()} ({self.channel_str}) {self.linked_object.s_identifier}\ndata integration time vs. onset distribution stats", 
                         fontsize=TITLE_FONTSIZE+9)
        else:
            ax.set_title(title, fontsize=TITLE_FONTSIZE+9)

        set_standard_ticks(ax=ax)
        # More tick settings (nonstandard)
        ax.tick_params(which="both", labelsize=TICK_LABELSIZE+9)
        hour_minute_format = DateFormatter("%H:%M")
        ax.yaxis.set_major_formatter(hour_minute_format)

        ax.grid(visible=grid, axis="both")

        set_legend(ax=ax, legend_loc="out", fontsize=LEGEND_SIZE+4)
        # ax.legend(loc=3, bbox_to_anchor=(1.0, 0.01), prop={'size': 24})

        if save:
            if not savepath:
                savepath = CURRENT_PATH
            plt.savefig(f"{savepath}{os.sep}int_time_vs_onset_distribution_stats_{self.spacecraft}_{self.sensor}_{self.linked_object.s_identifier}.png", transparent=False,
                        facecolor='white', bbox_inches='tight')

        plt.show()


    def show_onset_distribution(self, integration_time_index:int=0, xlim=None, show_background=True, 
                                save=False, savepath:str=None, legend_loc="out") -> None:
        """
        Displays all the unique onsets found with statistic_onset() -method in a single plot. The mode onset, that is
        the most common out of the distribution, will be shown as a solid line. All other unique onset times are shown
        in a dashed line, with the shade of their color indicating their relative frequency in the distribution.

        Note that statistic_onset() has to be run before this method. Otherwise KeyError will be raised.

        This is in all practicality just a dublicate method of Onset class' method with the same name. This exists just
        as an alternative to save computation time and use more memory instead.

        Parameters:
        -----------
        integration_time_index : int, default 0
                Choose which distribution from the integration time plot to show
        xlim : {tuple or list with __len__==2}, optional
                The limits of x-axis as pandas-compatible datetime strings, e.g., '2021-09-25 12:00'
        show_background : {bool}
                Boolean switch to show the background on the plot.
        save : {bool}, default False
                Boolean save switch.
        savepath : {str}, optional
                The directory path or subdirectory to save the figure.
        legend_loc : {str}, optional
                Either 'in' or 'out' to control wether the legend is inside or outside the figure borders.
        """

        rcParams["font.size"] = 20

        flux_series = self.list_of_series[integration_time_index]
        most_likely_onset = self.archive[integration_time_index]["most_likely_onset"]
        onsets = self.archive[integration_time_index]["unique_onsets"]

        # Gets the date of the event
        figdate = get_figdate(self.archive[0]["onset_list"])

        # Create a colormap for different shades of red for the onsets. Create the map with the amount of onsets+2, so that
        # the first and the final color are left unused. They are basically white and a very very dark red.
        cmap = plt.get_cmap("Reds", len(onsets)+2)

        # Creating the figure and axes
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        # x-axis settings:
        if not xlim:
            xlim = (most_likely_onset[0]-pd.Timedelta(hours=3), most_likely_onset[0]+pd.Timedelta(hours=2))
            ax.set_xlim(xlim)
        else:
            ax.set_xlim(pd.to_datetime(xlim[0]), pd.to_datetime(xlim[1]))

        set_standard_ticks(ax=ax)
        ax.xaxis_date()
        utc_dt_format1 = DateFormatter("%H:%M")
        ax.xaxis.set_major_formatter(utc_dt_format1)

        # This is what we draw and check y-axis limits against to
        flux_in_plot = flux_series.loc[(flux_series.index > xlim[0]) & (flux_series.index < xlim[1])]

        # y-axis: 
        ylim = set_fig_ylimits(ax=ax, ylim=None, flux_series=flux_in_plot)
        ax.set_ylim(ylim)

        ax.set_yscale("log")
        ax.set_ylabel(self.linked_object.unit, fontsize=AXLABEL_FONTSIZE)
        ax.set_xlabel(f"Time ({figdate})", fontsize=AXLABEL_FONTSIZE)

        # The intensity data, kw where dictates where the step happens. "mid" for the middle of the bin
        ax.step(flux_in_plot.index, flux_in_plot.values, where="mid")

        # Get the Onset object that was used to create this OnsetStatsArray, and call its BootstrapWindow to draw the background
        if show_background:
            self.linked_object.background.draw_background(ax=ax)

        # Marking the onsets as vertical lines of varying shades of red
        for i, onset in enumerate(onsets):

            # Not the most likely onset
            if onset[0] != most_likely_onset[0]:
                linestyle = (5, (10, 3)) if onset[1] > 0.2 else (0, (5, 5))
                ax.axvline(x=onset[0], c=cmap(i+1), ls=linestyle, label=f"{str(onset[0])[11:19]} ({np.round(onset[1]*100,2):.2f} %)")

            # The most likely onset
            # I'm accessing the i+1th color value in the map because I never want the first value (that is white in the standard Reds colormap)
            else:
                ax.axvline(x=onset[0], c="red", label=f"{str(onset[0])[11:19]} ({np.round(onset[1]*100,2):.2f} %)") 

        set_legend(ax=ax, legend_loc=legend_loc, fontsize=LEGEND_SIZE)

        int_time_str = f"{self.integration_times[integration_time_index]} integration time" if pd.Timedelta(self.integration_times[integration_time_index]) != self.linked_object.get_minimum_cadence() else f"{self.integration_times[integration_time_index]} data"
        # int_time_str = f"{self.integration_times[index]} integration time" if index != 0 else f"{int(self.linked_object.get_minimum_cadence().seconds/60)} min data" if self.linked_object.get_minimum_cadence().seconds>59 else f"{self.linked_object.get_minimum_cadence().seconds} s data"
        ax.set_title(f"{self.spacecraft.upper()}/{self.sensor.upper()} ({self.channel_str}) {self.linked_object.s_identifier}\nOnset distribution ({int_time_str})", fontsize=TITLE_FONTSIZE)

        if save:
            if not savepath:
                savepath = CURRENT_PATH
            plt.savefig(f"{savepath}{os.sep}onset_distribution_{self.spacecraft}_{self.sensor}_{self.linked_object.s_identifier}.png", transparent=False,
                        facecolor="white", bbox_inches="tight")

        plt.show()


    def show_onset_statistics(self, integration_time_index:int=0, xlim=None, show_background=True, 
                              save=False, savepath:str=None, legend_loc="out") -> None:
        """
        Shows the median, mode, mean and confidence intervals for the distribution of onsets got from statistic_onset().

        This is in all practicality just a dublicate method of Onset class' method with the same name. This exists just
        as an alternative to save computation time and use more memory instead.

        Parameters:
        -----------
        integration_time_index : int, default 0
                Choose which distribution from the integration time plot to show
        xlim : {tuple, list} of len()==2
                The left and right limits for the x-axis as pandas-compatible time strings.
        show_background : {bool}, optional
                Boolean switch to show the used background window on the plot
        save : {bool}, default False
                Boolean save switch.
        savepath : {str}, optional
                The directory path or subdirectory to save the figure to.
        """

        # Collect the median, mode and mean onsets from the distribution
        onset_median = self.archive[integration_time_index]["median_onset"]
        onset_mode = self.archive[integration_time_index]["most_likely_onset"][0]
        onset_mean = self.archive[integration_time_index]["mean_onset"]
        confidence_intervals1 = self.archive[integration_time_index]["1-sigma_confidence_interval"]
        confidence_intervals2 = self.archive[integration_time_index]["2-sigma_confidence_interval"]

        figdate = onset_mean.date().strftime("%Y-%m-%d")

        # This is just for plotting the time series with the chosen resolution
        flux_series = self.list_of_series[integration_time_index]

        # Plot commands and settings:
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        ax.step(flux_series.index, flux_series.values, where="mid")

        # Get the Onset object that was used to create this OnsetStatsArray, and call its BootstrapWindow to draw the background
        if show_background:
            self.linked_object.background.draw_background(ax=ax)

        # Vertical lines for the median, mode and mean of the distributions
        ax.axvline(onset_median, c=COLOR_SCHEME["median"], label="median")
        ax.axvline(onset_mode, c=COLOR_SCHEME["mode"], label="mode")
        ax.axvline(onset_mean, c=COLOR_SCHEME["mean"], label="mean")

        # 1-sigma uncertainty shading
        ax.axvspan(xmin=confidence_intervals1[0], xmax=confidence_intervals1[1], color=COLOR_SCHEME["1-sigma"], alpha=0.3, label="~68 % confidence")

        #2-sigma uncertainty shading
        ax.axvspan(xmin=confidence_intervals2[0], xmax=confidence_intervals2[1], color=COLOR_SCHEME["2-sigma"], alpha=0.3, label="~95 % confidence")


        # x-axis settings:
        if not xlim:
            xlim = (onset_median-pd.Timedelta(hours=3), onset_median+pd.Timedelta(hours=2))
            ax.set_xlim(xlim)
        else:
            ax.set_xlim(pd.to_datetime(xlim[0]), pd.to_datetime(xlim[1]))

        set_standard_ticks(ax=ax)
        ax.xaxis_date()
        utc_dt_format1 = DateFormatter("%H:%M")
        ax.xaxis.set_major_formatter(utc_dt_format1)

        # y-axis: 
        flux_in_plot = flux_series.loc[(flux_series.index > xlim[0]) & (flux_series.index < xlim[1])]
        ylim = set_fig_ylimits(ax=ax, ylim=None, flux_series=flux_in_plot)

        # If the lower y-boundary is some ridiculously small number, adjust the y-axis a little
        if np.log10(ylim[1])-np.log10(ylim[0]) > 10:
            ylim[0] = 1e-4

        ax.set_ylim(ylim)
        ax.set_yscale("log")
        ax.set_ylabel(self.linked_object.unit, fontsize=AXLABEL_FONTSIZE)
        ax.set_xlabel(f"Time ({figdate})", fontsize=AXLABEL_FONTSIZE)

        set_legend(ax=ax, legend_loc=legend_loc, fontsize=LEGEND_SIZE)

        # particle_str = "electrons" if self.species=='e' else "protons" if self.species=='p' else "ions"
        int_time_str = f"{self.integration_times[integration_time_index]} integration time" if pd.Timedelta(self.integration_times[integration_time_index]) != self.linked_object.get_minimum_cadence() else f"{self.integration_times[integration_time_index]} data"
        # int_time_str = f"{self.integration_times[index]} integration time" if index != 0 else f"{int(self.linked_object.get_minimum_cadence().seconds/60)} min data" if self.linked_object.get_minimum_cadence().seconds>59 else f"{self.linked_object.get_minimum_cadence().seconds} s data"
        ax.set_title(f"{self.spacecraft.upper()}/{self.sensor.upper()} ({self.channel_str}) {self.linked_object.s_identifier}\nOnset statistics ({int_time_str})", fontsize=TITLE_FONTSIZE)

        if save:
            if not savepath:
                savepath = CURRENT_PATH
            plt.savefig(f"{savepath}{os.sep}onset_statistics_{self.spacecraft}_{self.sensor}_{self.linked_object.s_identifier}.png", transparent=False,
                        facecolor='white', bbox_inches='tight')

        plt.show()


    def plot_cdf(self, index:int=0, prints:bool=False) -> None:
        """
        Plots the Cumulatice Distribution Function (CDF) for a distribution of onset times.

        This is in all practicality just a dublicate method of Onset class' method with the same name. This exists just
        as an alternative to save computation time and use more memory instead.

        Parameters:
        -----------
        index : int, default 0
                Choose which resampled distribution to plot
        prints : bool, default False
                If True, will also print out the median onset time and the confidence intervals.
        """

        # Remember that onset times are sorted according to their probability in ascending order
        onset_timestamps = np.array([pair[0] for pair in self.archive[index]["unique_onsets"]])
        probabilities = np.array([pair[1] for pair in self.archive[index]["unique_onsets"]])

        # Here let's sort the onsets and probabilities in temporal order:
        onset_time_indices = onset_timestamps.argsort() # numpy array.argsort() returns a list of indices of values in ascending order

        onset_timestamps = pd.to_datetime(onset_timestamps[onset_time_indices])
        probabilities = probabilities[onset_time_indices]

        # Effectively calculate the cumulative sum of probabilities
        cumulative_probabilities = [np.nansum(probabilities[:i+1]) for i in range(len(probabilities))]
        
        # Get the median, mode, mean and confidence intervals (1sigma and 2sigma of normal distribution) of the onset distribution
        median_onset, _, _, confidence_intervals_all = self.get_distribution_statistics(onset_statistics=self.archive[index], 
                                                                                    percentiles=[(15.89,84.10), (2.30, 97.7)])
        
        # Flatten the sigma lists for plotting. Look at this list comprehension :D
        confidence_intervals_all = [timestamp for sublist in confidence_intervals_all for timestamp in sublist]

        # Init the figure
        fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE)

        rcParams["font.size"] = 20

        int_time_str = f"{self.integration_times[index]} integration time" if pd.Timedelta(self.integration_times[index]) != self.linked_object.get_minimum_cadence() else f"{self.integration_times[index]} data"
        # int_time_str = f"{self.integration_times[index]} integration time" if index != 0 else f"{int(self.linked_object.get_minimum_cadence().seconds/60)} min data" if self.linked_object.get_minimum_cadence().seconds>59 else f"{self.linked_object.get_minimum_cadence().seconds} s data"

        # Settings for axes
        ax.xaxis_date()
        ax.set_xlabel("Time", fontsize=AXLABEL_FONTSIZE)
        ax.set_ylabel("Cumulative Probability", fontsize=AXLABEL_FONTSIZE)
        ax.set_title(f"Cumulative Distribution Function, ({int_time_str})", fontsize=TITLE_FONTSIZE)
        hour_minute_format = DateFormatter("%H:%M")
        ax.xaxis.set_major_formatter(hour_minute_format)

        # Plot command(s)
        ax.step(onset_timestamps, cumulative_probabilities, zorder=3, c='k', where="post")

        colors = ("red", "red", "blue", "blue")
        # Horizontal dashed lines for 1-sigma and 2-sigma limits
        for i, num in enumerate((0.159, 0.841, 0.023, 0.977)):
            label = "15.9-84.1 percentiles" if num==0.841 else ("2.3-97.7 percentiles" if num==0.977 else None)
            ax.axhline(y=num, ls="--", c=colors[i], zorder=1, alpha=0.8, label=label)

        # Vertical dashed lines for 1-sigma and 2-sigma limits
        for i, time in enumerate(confidence_intervals_all):
            ax.axvline(x=time, ls="--", c=colors[i], zorder=1, alpha=0.8)

        ax.axvline(x=median_onset, ls="--", c="maroon", zorder=2, alpha=0.8, label="median onset (50th percentile)")
        ax.axhline(y=0.50, ls="--", c="maroon", zorder=1, alpha=0.8)

        ax.legend(loc=10, bbox_to_anchor=(1.27, 0.88))

        if prints:
            print(f"The median onset time is: {str(median_onset.time())[:8]}")
            # print(f"The np.percentile()-calculated median onset time (50th percentile) is: {str(median_onset.time())[:8]}")
            display(Markdown("~68 % confidence interval"))
            print(f"{str(confidence_intervals_all[0].time())[:8]} - {str(confidence_intervals_all[1].time())[:8]}")
            print(" ")
            display(Markdown("~95 % confidence interval"))
            print(f"{str(confidence_intervals_all[2].time())[:8]} - {str(confidence_intervals_all[3].time())[:8]}")

        plt.show()


    def calculate_weighted_uncertainty(self, weight_type:str="inverse_variance", returns=False):
        """
        Calculates the weighted confidence intervals based on the confidence intervals of the varying
        integration times. Weighting is done so that the confidence intervals are in a way normalized to 
        to the confidence intervals acquired from the native data resolution.

        Parameters:
        -----------
        weight_type : {str} either 'int_time', 'uncertainty' or 'inverse_variance' (default)
                    Defines the logic at which uncertainties are weighted.

        returns : {bool} default False
                    Switch for this method to also return the weights.
        """

        if weight_type not in ("uncertainty", "int_time", "inverse_variance"):
            raise ValueError(f"Argument {weight_type} is not a valid value for weight_type! It has to be either 'uncertainty', 'int_time' or 'inverse_variance'.")

        # Asserting the weights so that w_0 = integration time of final distribution and w_{-1} = 1. 
        int_weights = [w for w in range(len(self.archive), 0, -1)]

        # Collect the confidence intervals and median to their respective arrays
        sigma1_low_bounds = np.array([stats["1-sigma_confidence_interval"][0] for stats in self.archive])
        sigma2_low_bounds = np.array([stats["2-sigma_confidence_interval"][0] for stats in self.archive])
        modes = np.array([stats["most_likely_onset"][0] for stats in self.archive])
        medians = np.array([stats["median_onset"] for stats in self.archive])
        sigma1_upper_bounds = np.array([stats["1-sigma_confidence_interval"][1] for stats in self.archive])
        sigma2_upper_bounds = np.array([stats["2-sigma_confidence_interval"][1] for stats in self.archive])

        if weight_type == "inverse_variance":

            # Inverse-variance weighting:
            # https://en.wikipedia.org/wiki/Inverse-variance_weighting
            weights = self.inverse_variance_weights()

        elif weight_type == "uncertainty":

            # Instead weighting by the width of 2-sigma uncertainty intervals
            sigma2_widths = np.array([(sigma2_upper_bounds[i] - low_bound) for i, low_bound in enumerate(sigma2_low_bounds)])

            # Convert the widths to minutes as floating point numbers
            sigma2_widths = np.array([width.seconds/60 for width in sigma2_widths])

            # The weights here are the inverse's of the widths of the 2-sigma uncertainty intervals
            weights = np.array([1/width for width in sigma2_widths])

            if len(weights)==0:
                weights = np.array(int_weights)

        # Using integration time as weighting
        else:
            weights = np.array(int_weights)

        # Weighting the confidence intervals and the median
        self.w_sigma1_low_bound = weight_timestamp(weights=weights, timestamps=sigma1_low_bounds)
        self.w_sigma2_low_bound = weight_timestamp(weights=weights, timestamps=sigma2_low_bounds)
        self.w_mode = weight_timestamp(weights=weights, timestamps=modes)
        self.w_median = weight_timestamp(weights=weights, timestamps=medians)
        self.w_sigma1_upper_bound = weight_timestamp(weights=weights, timestamps=sigma1_upper_bounds)
        self.w_sigma2_upper_bound = weight_timestamp(weights=weights, timestamps=sigma2_upper_bounds)

        # Introduce a check to move the boundaries away from mode and median if they are too close
        self.check_weighted_timestamps()

        # Propagate the infromation of the weighted median and confidence intervals to the linked object too
        self.set_w_median_and_confidence_intervals( self.w_mode,
                                                    self.w_median,
                                                    self.w_sigma1_low_bound, self.w_sigma1_upper_bound,
                                                    self.w_sigma2_low_bound, self.w_sigma2_upper_bound)

        if returns:
            return self.w_mode, self.w_median, self.w_sigma1_low_bound, self.w_sigma1_upper_bound, self.w_sigma2_low_bound, self.w_sigma2_upper_bound

    def check_weighted_timestamps(self):
        """
        Checks that the intervals boundaries are separated from the mode and the median by at least half
        of the finest cadence used (minimum integration time). Moves them if not.
        """

        # We will not allow for a separation of less than the smalles integration time used
        #min_separation = self.linked_object.get_minimum_cadence()/2
        min_separation = pd.Timedelta(self.integration_times[0])/2

        if self.w_mode - self.w_sigma1_low_bound < min_separation:
            self.w_sigma1_low_bound = self.w_mode - min_separation
        if self.w_mode - self.w_sigma2_low_bound < min_separation:
            self.w_sigma2_low_bound = self.w_mode - min_separation
        
        if self.w_sigma1_upper_bound - self.w_mode < min_separation:
            self.w_sigma1_upper_bound = self.w_mode + min_separation
        if self.w_sigma2_upper_bound - self.w_mode < min_separation:
            self.w_sigma2_upper_bound = self.w_mode + min_separation


    def inverse_variance_weights(self):
        """
        Calculates the inverse variance weights for onset distributions.
        
        Uses the inverse variance IN MINUTES as weights.
        
        Returns: weights : {np.ndarray} array of shape==(len(archive),)
        """

        # This is the variance of a distribution with 1000 onset times wide 1 minute equally spaced.
        # Experimentally found with: pd.date_range("2024-01-01 00:00:00", "2024-01-01 00:01:00", 1000)
        VAR_OF_MINUTE_WIDE_DISTRIBUTION = 0.08350016683180035

        weights = np.zeros([len(self.archive)], dtype=float)

        # Loop through each archive individually
        for i, stats in enumerate(self.archive):

            # Fetch onset list
            onsets = pd.DatetimeIndex(stats["onset_list"])

            # Convert to nanoseconds, because variance can NOT be calculated for
            # datetime type objects
            onsets_as_nanoseconds = onsets.astype(np.int64)

            # NaTs equal to negative nanosecond values. Replace with nans as they are 
            # not considered for calculation anyway
            onsets_as_nanoseconds = np.where(onsets_as_nanoseconds<0, np.nan, onsets_as_nanoseconds)

            # Calculating the variance of these numbers is in essence the variance
            # of the onset times (in units of nanoseconds squared)
            variance_in_nanoseconds = np.nanvar(onsets_as_nanoseconds)

            # Multiply with 10^-18 to convert nanoseconds squared to seconds squared
            variance_in_seconds = variance_in_nanoseconds * 1e-18

            # Dividing by 3600 = 60*60 converts seconds squared to minutes squared
            variance_in_minutes = variance_in_seconds / 3600.

            # This check is here because sometimes numpy calculates the variance
            # for a delta function to be an extremely small number of the order of
            # 10^-17. I suspect it is due to calculating precision as in the case of
            # a delta function the variance is 0. Anyway replace that here with 
            # the variance of an equi-spaced distribution with width 1 minute times the time reso.
            weights[i] = variance_in_minutes if variance_in_minutes > 1e-10 else VAR_OF_MINUTE_WIDE_DISTRIBUTION * (i+1)

        # Return the reciprocal (1 / weights) 
        return np.reciprocal(weights)



def check_confidence_intervals(confidence_intervals, time_reso):
    """
    Small helper function to make sure that the confidence intervals are at least as large as the time resolution
    of the time series data.
    """

    # Initialize minimum inner boundaries (no outer boundary can be inside these)
    min_low_boundary = confidence_intervals[0][0]
    min_high_boundary = confidence_intervals[0][1]

    new_intervals = []
    for i, interval in enumerate(confidence_intervals):

        if interval[1] - interval[0] < pd.Timedelta(time_reso):

            # Calculate the center of this interval
            center_of_interval = datetime_mean(np.array([interval[0], interval[1]]))

            # Push the boundaries of the interval further from their center by half of the time resolution
            new_interval0, new_interval1 = center_of_interval - pd.Timedelta(time_reso)/2, center_of_interval + pd.Timedelta(time_reso)/2

            # Check that inner boundaries are never outside outer boundaries
            # Also check that 1sigma boundary is never outside 2sigma boundaries
            if i==0:
                min_low_boundary = new_interval0
                min_high_boundary = new_interval1

            if new_interval0 > min_low_boundary:
                new_interval0 = min_low_boundary
            
            if new_interval1 < min_high_boundary:
                new_interval1 = min_high_boundary

        # Even if the interval length is larger than the time resolution, some outer boundary may be left inside
        # an inner interval
        else:

            # Here check that outer interval boundaries are not inside inner interval boundaries
            if interval[0] > min_low_boundary:
                new_interval0 = min_low_boundary
            else:
                new_interval0 = interval[0]
            
            if interval[1] < min_high_boundary:
                new_interval1 = min_high_boundary
            else:
                new_interval1 = interval[1]

        new_intervals.append((new_interval0,new_interval1))

    return new_intervals


def weight_timestamp(weights, timestamps):
    """
    Calculates weighted timestamp from a list of weights and timestamps.

    Parameters:
    -----------
    weights : {list of weights}
    timestamp : {list of pd.datetimes}
    """

    # It could happen that there is only one timestamp and that it's a NaT. In this case numpy will falsely
    # average the NaT to 1970-01-01 00:00, so check it here to avert this.
    if len(timestamps)==1 and isinstance(timestamps[0], pd._libs.tslibs.nattype.NaTType):
        return pd.NaT

    # First make sure timestamps come numpy datetime64 format (in nanoseconds)
    timestamps = np.array([t.asm8 for t in timestamps])

    # Mask the weights and timestamps that contain nans/NaTs with these indices
    mask = ~np.isnan(timestamps)

    # Calculate weighted_avg. Use np.ndarray.view() to convert datetime values to floating point numbers so that
    # math can be done with them. The timestamp that this yields will also be a floating point number
    wavg_ns = np.average(timestamps[mask].view(dtype="float64"), weights=weights[mask])

    # Convert the float to a numpy timestamp
    wavg_timestamp = np.array(wavg_ns).view(dtype="datetime64[ns]")

    # Return the weighted average timestamp (final conversion to pandas datetime is perhaps not necessary, but let's be sure)
    return pd.to_datetime(wavg_timestamp)


def subwindows(Window, window_placement="equal", fraction:float = 0.2):
    """
    A function that turns a BootstrapWindow object into a list of smaller BootstrapWindows according to n_shifts.
    Each new smaller BootstrapWindow will still be contained within the original Window.
    Returned windows fill always have n_shifts set to 0.

    Parameters:
    -----------
    Window : BootstrapWindow object
            Defines the boundaries within the small new windows will be created, and the amount of them
    window_placement : str
            Defines the type of subwindow distribution. Either 'random' or 'equal'
    fraction : float
            The size fraction of the input window that the subwindows will be. Defaults to 0.2, i.e., 20 %
    Returns:
    -----------
    windowlist : list(BootstrapWindow)
            A list of <n_shifts> new Bootstrapwindows.
    """

    # This is the amount of samples that are picked from the big window, use the same for the little windows
    bootstraps = Window.bootstraps

    # We need to know that temporal range of the big window, because the little windows have to be some reasonable
    # fraction of the big window
    big_windowrange_in_secs = (Window.end - Window.start).seconds

    if window_placement == "random":

        # Small window will reach for fraction times the big window's range
        small_window_len = int(fraction * big_windowrange_in_secs)

        # Pre-calculate random starting points for the little windows, remembering that they cannot start later
        # than (1 - fraction) of the big window, in order for the big window to completely contain all the little windows
        window_startpoints = (1.0 - fraction) * big_windowrange_in_secs * np.random.random((Window.n_shifts,))

    elif window_placement == "equal":

        # Small window will in the equal case be 1/n_shifts the size of the big windows, so as to fit equal-sized
        # windows inside the big window
        fraction = 1/Window.n_shifts
        small_window_len = int(fraction * big_windowrange_in_secs)

        # Generate starting points for the windows, we already know everything else about them
        window_startpoints = [sec for sec in np.linspace(start=0, stop=big_windowrange_in_secs-small_window_len, num=Window.n_shifts)]

    else:
        raise Exception(f"The window_placement must be either 'random' or 'equal'. The input was: {window_placement}")


    windowlist = []
    for sec in window_startpoints:

        start = Window.start + pd.Timedelta(seconds=int(sec))
        end = start + pd.Timedelta(seconds=small_window_len)

        # New bootstrapwindows is initialized with random start, end that's at a set distance from the start,
        # n_shifts always 0 and the amount of bootstrapped samples is equal to the original amount
        new_window = BootstrapWindow(start=start, end=end, n_shifts=0, bootstraps=bootstraps)

        windowlist.append(new_window)

    return windowlist


def get_x_errors(e_min, e_max, inverse_betas, mass_energy:float):
        """
        Calculates the width of energy channels in terms of the inverse speed (beta = v/c)

        Parameters:
        -----------
        e_min : array-like
                    The lower boundaries of energy channels in terms of eVs
        e_max : array-like
                    The higher boundaries of energy channels in terms of eVs
        inverse_betas : array-like
                    The inverse betas corresponding to channel nominal energies
        mass_energy : float
                    E = mc^2, the energy related to the mass of the particle

        Returns:
        -----------
        x_errors_lower : tuple
                    The lower boundaries of energy channels in terms of inverse beta
        x_errors_upper : tuple
                    The higher boundaries of energy channels in terms of inverse beta
        x_errors : tuple
                    The widths of energy channels in terms of 1/beta
        """

        # Error bars in x direction:
        # First turn channel energy boundaries to Joules
        e_min_joule = [((e*u.eV).to(u.J)).value for e in e_min]
        e_max_joule = [((e*u.eV).to(u.J)).value for e in e_max]

        # Calculate beta from channel energy boundaries
        beta_lower = [np.sqrt(1-(( e_j/mass_energy + 1)**(-2))) for e_j in e_min_joule]
        beta_upper = [np.sqrt(1-(( e_j/mass_energy + 1)**(-2))) for e_j in e_max_joule]

        # The Inverse of beta
        beta_inv_lower = [1/x for x in beta_lower]
        beta_inv_upper = [1/x for x in beta_upper]

        # Calculate the difference between boundary and actual value
        x_errors_lower = [x-inverse_betas[indx] for indx, x in enumerate(beta_inv_lower)]
        x_errors_upper = [inverse_betas[indx]-x for indx, x in enumerate(beta_inv_upper)]

        x_errors = [x-x_errors_upper[i] for i, x in enumerate(x_errors_lower)]

        return x_errors_lower, x_errors_upper, x_errors


def get_time_errors(onset_times, spacecraft):
    """
    Gets the errors in y-direction (time).

    Parameters:
    -----------
    onset_times : array-like
                All the onset times
    spacecraft : str
                The spacecraft identifier string.
    
    Returns:
    -----------
    time_error : list(pandas.Timedelta)
                A list of time errors as pandas-Timedelta objects
    """

    # Error in the y-direction (time):
    time_error = []

    for ts in onset_times:

        try:

            if isinstance(ts.freqstr,str):

                # freqstr could be negative, check!
                if ts.freqstr[0] == '-':
                    freqstr = ts.freqstr[1:]
                else:
                    freqstr = ts.freqstr
                time_error.append(pd.Timedelta(freqstr))

            else:
                if spacecraft in ["stereo-a", "stereo a", "stereo-b", "stereo b", "sta", "stb"]:
                    time_error.append(pd.Timedelta('1T'))
                else:
                    time_error.append(pd.Timedelta('1S'))

        except AttributeError:

            if spacecraft in ["stereo-a", "stereo a", "stereo-b", "stereo b", "sta", "stb"]:
                time_error.append(pd.Timedelta('1T'))
            else:
                time_error.append(pd.Timedelta('1S'))

    return time_error


def calculate_inverse_betas(channel_energies, mass_energy):
    """
    Calculates the inverse unitless speeds (1/beta = c/v) of particles given channel energies.

    Parameters:
    -----------
    channel_energies : array-like
                The nominal energies of energy channels.
    
    Returns:
    -----------
    inverse_betas : list
                The inverse betas corresponding to the nominal energies of the channels.
    """

    # Transform channel energies in eV to Joules
    e_Joule = [((En*u.eV).to(u.J)).value for En in channel_energies]

    # Beta, the unitless speed (v/c)
    beta = [np.sqrt(1-(( e_j/mass_energy + 1)**(-2))) for e_j in e_Joule]

    # Invert the beta values
    inverse_betas = np.array([1/b for b in beta])

    return inverse_betas


def tsa(t0:datetime.datetime, L:float, v:float):
    """
    Calculates the time shifted injection time assuming a path length L and particle speed v

    Parameters:
    -----------
    t0 : datetime.datetime
            The onset time of the event
    L : float
            The path length in AU
    v : float
            particle speed in m/s
    """

    path_in_meters = L * u.au.to(u.m)
    v_in_ms = v * u.m / u.s
    
    path_delta_sec = path_in_meters / v_in_ms

    t_inj = t0 - pd.Timedelta(seconds=path_delta_sec.value)

    return t_inj


def datetime_to_sec(onset_times):
    """
    Returns the the timestamps of onset times as seconds since the linux Epoch (first of January, 1970)

    Parameters:
    -----------
    onset_times : {array_like}
    """

    # date_in_sec is the list that goes to the fitting function
    dates_in_sec = []

    # Calculate the difference between an onset time and the start of the epoch in seconds
    for date in onset_times:
        dates_in_sec.append((date.to_pydatetime() - datetime.datetime(1970,1,1,0,0,0)).total_seconds())

    return dates_in_sec


def datetime_nanmedian(timestamps):
    """
    Finds the median datetime, ignoring NaTs.

    timestamps : {array-like}
                    An array or list of pandas-compatible timestamps.
    """

    if not isinstance(timestamps,np.ndarray):
        timestamps = np.array(timestamps)

    # Mask away NaTs. pd.isnull() also applies on NaTs
    timestamps_masked = timestamps[~pd.isnull(timestamps)]

    # Sort the array
    timestamps_masked.sort()
    
    # Pick the middlemost value to return
    return timestamps_masked[len(timestamps_masked)//2]


def sample_mean_and_std(start, end, flux_series, sample_size=None, prints_warning=True):
    """
    Calculates the mean and standard deviation of the background period by
    taking a random sample of the background window of size sample_size * background_window.

    Parameters:
    -----------
    start: datetime
            Starting point of the calculation
    end: datetime
            Ending point of the calculation
    flux_series: pandas Series
            An array of observed values, indexed by datetime
    sample_size: float, default None
            A fraction of all the points to be considered in the calculation
            of mean and std. The sample is a random sample.
    prints_warning: bool, default True
            The function prints a warning if the sample size is small (data points < 100)

    Returns:
    --------
    (mean, std): tuple
            mean is the mean of the background period, std is its standard deviation
    """

    # Background contains all the points that could be considered for the calculation
    # of mean and std
    background = flux_series.loc[(flux_series.index >= start) & (flux_series.index < end)]

    # If sample_size is not provided, assume that all data points are picked
    if not sample_size:
        print("sample_size was not provided, assuming a 100 % sample.")
        sample_size = int(len(background))

    # Sample size is a fraction of the amount of background window data points
    else:
        if sample_size < 0.0 or sample_size > 1.0:
            raise Exception(f"The sample size must be between 0 and 1, not {sample_size}!")
        else:
            sample_size = int(sample_size * len(background))
            if len(background) < 100 and prints_warning:
                print(f"Warning, random sample is only {sample_size} data points!")


    # A random sample of all the background data points.
    sample = np.random.choice(background, replace=False, size=sample_size)

    # Calculate mean and std of the chosen sample of data points
    mean = np.nanmean(sample)
    std = np.nanstd(sample)

    return (mean, std)


def detrend_onsets(timestamps, shift):
    """
    Takes as an input an array of datetimes, and shifts them forward in time by the given amount.
    
    Parameters:
    ---------
    timestamps : array-like
                A list or an array containing the timestamps to be shifted
    shift : str
                Pandas-compatible timedelta string to shift the timestamps, e.g., '30s' or '1min'
    
    Returns:
    --------
    new_timestamps : array-like
                The detrended timestamps
    """

    # Init a new array for the timestamps to ensure that vectorization works
    timestamps_arr = np.asarray(timestamps, dtype=object)
    new_timestamps = timestamps_arr + pd.Timedelta(shift)

    return new_timestamps


def datetime_mean(arr):
    """
    Returns the mean of an array of datetime values
    """
            
    arr1 = pd.Series(arr)
    
    return arr1.mean()


def find_biggest_nonzero_unit(timedelta):
    """
    Finds the biggest unit of time that describes the timedelta.
    Available units of time include a day, an hour, a minutes and a second.
    """
    
    days, hours, minutes, seconds, _, _, _ = timedelta.components
    
    if days != 0:
        return 'D'
    if hours != 0:
        return 'H'
    if minutes != 0:
        return "min"
    return 's'


def find_earliest_possible_onset(flux_series, cusum_function, onset_time):
    """
    Finds the earliest possible onset time by backtracking the cusum_function starting from the onset time and finding the
    first nonzero occurence of the cusum function.

    Parameters:
    -----------
    flux_series : pandas series
                    The series from which the onset was found. Intensity measurements indexed by datetimes.
    cusum_function : array-like
                    The values of the cusum function that is used to identify the onset time in the aforementioned time series.
    onset_time : datetime
                    The true onset time found from the time series with the cusum method.
    Returns:
    -----------
    earliest_possible_onset : datetime
                    The earliest possible onset time, identified by the first nonzero occurence of the cusum function 
                    backtracking from the true onset time.
    """
    # First to find the index of the onset, use numpy.where()
    datetimes = flux_series.index

    index_of_onset = np.where(datetimes == onset_time)[0][0]

    # Before looping througfh the cusum_function, flip it to make sure we are advancing in the past, not towards the future
    flipped_cusum_function = np.flip(cusum_function)

    # Init the earliest possible onset time
    earliest_possible_onset = onset_time

    # We'll go though the cusum_function, backtracking towards the past starting from the onset time
    for i in range(index_of_onset, len(flipped_cusum_function)):

        # The earliest possible onset time is the first nonzero occurence of cusum function PRIOR to the true onset
        if flipped_cusum_function[i] == 0:
            earliest_possible_onset = datetimes[i-1]
            break

    return earliest_possible_onset




# ==============================================================================================

def get_time_reso(series) -> str:
    """
    Returns the time resolution of the input series.

    Parameters:
    -----------
    series: Pandas Series

    Returns:
    ----------
    resolution: {str}
            Pandas-compatible freqstr
    """

    if series.index.freq is None:
        
        # Better approach:
        index_diffs = series.index.diff()
        
        # There might be unregular time differences, pick the most
        # appearing one -> mode.
        diffs, counts = np.unique(index_diffs, return_counts=True)
        mode_dt = pd.Timedelta(diffs[np.argmax(counts)])

        # Round up to the nearest second, because otherwise e.g., STEREO / SEPT data
        # that may have cadence of '59.961614005' seconds is interpreted to have nanosecond
        # precision.
        mode_dt = mode_dt.round(freq='s')

        # If less than 1 minute, express in seconds
        divisor = 60 if mode_dt.resolution_string == "min" else 1
        
        return f"{mode_dt.seconds//divisor} {mode_dt.resolution_string}"

    
    else:
        freq_str = series.index.freq.freqstr
        return freq_str if freq_str!="min" else f"1 {freq_str}"

# ===========================================================================================

def calc_chnl_nominal_energy(e_min, e_max, mode='gmean'):
    """
    Calculates the nominal energies of each channel, given the lower and upper bound of
    the channel energies. Channel energy should be given in eVs.

    Parameters:
    -----------
    e_min, e_max: np.arrays
            Contains lower and higher energy bounds for each channel respectively in eV
    mode: str, default='gmean'
        The mode for calculating nominal channel energy. Use 'gmean' for geometric mean, and
        'upper' for the upper energy bound.
    
    Returns:
    --------
    energies: np.array 
            Contains the nominal channel energy for each channel
    """

    # geometric mean
    if mode=='gmean':
        energies = np.sqrt(np.multiply(e_min,e_max))
        return energies
    
    # upper bound of energy channels
    elif mode=='upper':
        energies = e_max
        return energies
    
    else:
        raise Exception("Choose either 'gmean' or 'upper' as the nominal energy.")

#===========================================================================================

def path_length_calculator(distance, solar_wind = 400):
    """
    Calculates the length of the nominal parker spiral arm.

    Parameters:
    ------------
    distance : float, int
                heliospheric distance of the object in AU
    solar_wind = float, int
                solar wind speed in km/s, default = 400

    Returns:
    -----------
    path_length : float
                the path length in AU assuming the nominal Parker spiral
    """
    from sunpy.sun import constants as sunc

    # init the relevant constants
    OMEGA = sunc.sidereal_rotation_rate.value * (u.deg/u.day).to(u.rad/u.s) #rad/s, the equatorial sidereal rotation period of the Sun

    solar_wind_au = solar_wind * u.km.to(u.au)

    # These constants appears many times in the equation
    a = solar_wind_au/OMEGA
    psi = distance/a

    sqrt_term = np.sqrt(psi*psi + 1)

    # Calculate the length of the spiral arc
    path_length = 0.5 * a * ( psi*sqrt_term + np.log(psi + sqrt_term) ) # np.log() is the natural logarithm

    return path_length

#===========================================================================================

def flux2series(flux, dates, cadence=None):
    """
    Converts an array of observed particle flux + timestamps into a pandas series
    with the desired cadence.

    Parameters:
    -----------
    flux: an array of observed particle fluxes
    dates: an array of corresponding dates/times
    cadence: str - desired spacing between the series elements e.g. '1s' or '5min'

    Returns:
    ----------
    flux_series: Pandas Series object indexed by the resampled cadence
    """

    # from pandas.tseries.frequencies import to_offset

    # set up the series object
    flux_series = pd.Series(flux, index=dates)
    
    # if no cadence given, then just return the series with the original
    # time resolution
    if cadence is not None:
        try:
            flux_series = flux_series.resample(cadence, origin='start').mean()
            flux_series.index = flux_series.index + pd.tseries.frequencies.to_offset(pd.Timedelta(cadence)/2)
        except ValueError:
            raise Warning(f"Your 'resample' option of [{cadence}] doesn't seem to be a proper Pandas frequency!")

    return flux_series

# ======================================================================

def calculate_mean_and_std(start, end, flux_series):
    """
    This function calculates the mean and standard deviation of the background period
    which is used in the onset analysis. 

    Parameters:
    -----------
    start: datetime
            Starting point of the calculation
    end: datetime
            Ending point of the calculation
    flux_series: pandas Series
            An array of observed values, indexed by datetime

    Returns:
    --------
    (mean, std): tuple
            mean is the mean of the background period, std is its standard deviation
    """

    background = flux_series.loc[(flux_series.index >= start) & (flux_series.index < end)]

    mean = np.nanmean(background)
    std = np.nanstd(background)

    return (mean, std)

# =======================================================================

def erase_glitches(series, glitch_threshold, time_barr):
    """
    Every value above glitch_threshold before time_barr will be replace with nan.
    Also return the erased values in a series called "glitches", in case we need them for something, e.g. plotting.

    Parameters:
    -------------
    series : Pandas Series
                A series of intensity values indexed by pandas-compatible datetime-objects.
    glitch_threshold : float
                The maximum number above which in the right timeframe the values of the series are considered glitches.
    time_barr : datetime
                The time at which glitch_threshold stop applying.
    
    Returns:
    -------------
    series2 : Pandas Series
                The new series, with glitches omitted.
    glitches : Pandas Series
                All the found glitches collected to their own series, indexed by datetime.
    """

    if isinstance(time_barr,str):
        time_barr = pd.to_datetime(time_barr)

    # logically: if a value is BOTH before time_barr AND larger than glitch_threshold,
    # then the statement is false, and the corresponding value will be replaced with nan.
    glitches = series.loc[(series.values >= glitch_threshold) & (series.index < time_barr)]
    series2 = series.where((series.values < glitch_threshold) | (series.index >= time_barr))

    return series2, glitches

# ============================================================================

def onset_determination(ma_sigma, flux_series, cusum_window, avg_end, sigma_multiplier : int = 2) -> list :
    """
    Calculates the CUSUM function to find an onset time from the given time series data.

    Parameters:
    -----------
    ma_sigma : tuple(float, float)

    flux_series : pandas Series

    cusum_window : int

    avg_end : pandas datetime

    sigma_multiplier : float, int

    Returns:
    --------
    list(ma, md, k_round, h, norm_channel, cusum, onset_time)
    """

    # assert date and the starting index of the averaging process
    date = flux_series.index

    # First cut the part of the series that's before the ending of the averaging window. Then subtract that
    # from the size of the original series to get the numeral index that corresponds to the ending of the 
    # averaging window. starting_index is then the numeral index of the first datapoint that belongs outside the averaging window
    cut_series = flux_series.loc[flux_series.index > avg_end]
    start_index = flux_series.size - cut_series.size

    ma = ma_sigma[0]
    sigma = ma_sigma[1]

    md = ma + sigma_multiplier*sigma

    # k may get really big if sigma is large in comparison to mean
    try:
        k = (md-ma)/(np.log(md)-np.log(ma))

        # If ma == 0, then std == 0. Hence CUSUM should not be restricted at all -> k_round = 0
        # Otherwise k_round should be 1
        if not np.isnan(k):
            k_round = round(k/sigma) if k/sigma > 1 else k/sigma
        else:
            k_round = 1 if ma > 0 else 0

    except (ValueError, OverflowError) as error:
        # the first ValueError I encountered was due to ma=md=2.0 -> k = "0/0"
        # OverflowError is due to k = inf
        # print(error)
        k_round = 1 if ma > 0 else 0

    # choose h, the variable dictating the "hastiness" of onset alert
    h = 2 if k_round>1 else 1

    alert = 0
    cusum = np.zeros(len(flux_series))
    norm_channel = np.zeros(len(flux_series))
    
    # set the onset as default to be NaT (Not a daTe)
    onset_time = pd.NaT

    # start at the index where averaging window ends
    for i in range(start_index,len(cusum)):

        # normalize the observed flux
        norm_channel[i] = (flux_series[i]-ma)/sigma

        # calculate the value for ith cusum entry
        cusum[i] = max(0, norm_channel[i] - k_round + cusum[i-1])

        # check if cusum[i] is above threshold h, if it is -> increment alert
        if cusum[i]>h:
            alert=alert+1
        else:
            alert=0

        # cusum_window(default:30) subsequent increments to alert means that the onset was found
        if alert == cusum_window:
            onset_time = date[i - alert]
            break

    # ma = mu_a = background average
    # md = mu_d = background average + 2*sigma
    # k_round = integer value of k, that is the reference value to poisson cumulative sum
    # h = 1 or 2,describes the hastiness of onset alert
    # onset_time = the time of the onset

    return [ma, md, k_round, h, norm_channel, cusum, onset_time]


def onset_determination_cr(ma_sigma, flux_series, cusum_window, avg_end, sigma_multiplier : int = 2) -> list :
    """
    Calculates the CUSUM function to find an onset time from a given count rate data.
    
    The essential difference to the 'onset_determination()' -function is that this does NOT
    employ z-standardized intensity. rather it simply uses the counts for calculating
    CUSUM. Also k is not normalized to standard deviation. 

    Parameters:
    -----------
    ma_sigma : tuple(float, float)

    flux_series : pandas Series

    cusum_window : int

    avg_end : pandas datetime

    sigma_multiplier : float, int

    Returns:
    --------
    list(ma, md, k_round, h, flux_series.values, cusum, onset_time)
    """

    # Assert date and the starting index of the averaging process
    date = flux_series.index

    # First cut the part of the series that's before the ending of the averaging window. Then subtract that
    # from the size of the original series to get the numeral index that corresponds to the ending of the 
    # averaging window. starting_index is then the numeral index of the first datapoint that belongs outside
    # the averaging window
    cut_series = flux_series.loc[flux_series.index > avg_end]
    start_index = flux_series.size - cut_series.size

    ma = ma_sigma[0]
    sigma = ma_sigma[1]

    md = ma + sigma_multiplier*sigma

    # k may get really big if sigma is large in comparison to mean
    try:
        k = (md-ma)/(np.log(md)-np.log(ma))

        # If ma == 0, then std == 0. Hence CUSUM should not be restricted at all -> k_round = 0
        # Otherwise k_round should be 1
        if not np.isnan(k):
            k_round = round(k) if k > 1 else k
        else:
            k_round = 1 if ma > 0 else 0

    except (ValueError, OverflowError) as error:
        # the first ValueError I encountered was due to ma=md=2.0 -> k = "0/0"
        # OverflowError is due to k = inf
        # print(error)
        k_round = 1 if ma > 0 else 0

    # Choose h, the variable dictating the "hastiness" of onset alert
    h = 2 if k_round>1 else 1

    alert = 0
    cusum = np.zeros(len(flux_series))
    
    # Set the onset as default to be NaT (Not a daTe)
    onset_time = pd.NaT

    # Start at the index where averaging window ends
    for i in range(start_index,len(cusum)):

        # Calculate the value for ith cusum entry
        cusum[i] = max(0, flux_series[i] - k_round + cusum[i-1])

        # check if cusum[i] is above threshold h, if it is -> increment alert
        if cusum[i]>h:
            alert += 1
        else:
            alert = 0

        # cusum_window(default:30) subsequent increments to alert means that the onset was found
        if alert == cusum_window:
            onset_time = date[i - alert]
            break

    # ma = mu_a = background average
    # md = mu_d = background average + 2*sigma
    # k_round = integer value of k, that is the reference value to poisson cumulative sum
    # h = 1 or 2,describes the hastiness of onset alert
    # onset_time = the time of the onset

    return [ma, md, k_round, h, flux_series.values, cusum, onset_time]


#==============================================================================================

def bootstrap_mus_and_sigmas(flux_series, window_start, window_end, n_bstraps, sample_size=None):
    """
    The function looks at flux_series through a time window, and extracts n_bstraps samples
    from it. The size of a sample size <sample_size> must not be greater than the amount of data points
    inside the time window. It then calculates the mean and standard deviation of those samples and returns
    an array of means and an array of stds.
    """ 
    
    # make a selection including only the measurements in the window
    window = flux_series.loc[(flux_series.index>=window_start) & (flux_series.index<window_end)]

    if sample_size is None:
        #default sample size is set to 50 % of total datapoints in a window
        sample_size = int(len(window)/2)
    
        # Print out a warning if the amount of datapoints in a window is less than 100
        # This is done ONLY IF no sample size is provided
        if len(window)<100:
            print(f"Warning! The ensemble of background values is {len(window)}<100. \nRaising sample size to 75 %.")
            sample_size = int(3*len(window)/4)

    # if sample_size is not None, check that it's a reasonable value
    else:
        if sample_size < 0.0 or sample_size > 1.0:
            raise Exception(f"The sample size must be between 0 and 1, not {sample_size}!")
        else:
            sample_size = int(sample_size * len(window))
            if len(window) < 100:
                print(f"Warning, random sample is only {sample_size} data points!")

    # initialize arrays
    mean_arr = np.zeros(n_bstraps)
    std_arr = np.zeros(n_bstraps)

    # repeat n_bstraps times
    for i in range(n_bstraps):

        # look inside the window and randomly pick sample_size amount of measurements there
        sample = np.random.choice(window, replace=True, size=sample_size)
        
        # then calculate the mean and standard deviation of that sample
        mean_arr[i] = np.nanmean(sample)
        std_arr[i] = np.nanstd(sample)


    return mean_arr, std_arr

#================================================================================

def unique_vals(onset_list, prints=False):
    '''
    Returns all the unique values from a list of onset times, and the corresponding counts. Counts also include
    the amount of NaT values at the end of the array, which is why len(uniques) == len(counts)-1.

    Parameters:
    -----------
    onset_list : list
                A list of onset times
    prints : bool, default False
                Switch to print out the relative abundances of each onset in the list
    
    Returns:
    -----------
    uniques : np.array
                An array of the UNIQUE appearances of the onset times in the input list
    counts : np.array
                An array of the number of appearances of each onset in the input list + the amount of NaTs.
    '''
    
    # Convert to pandas datetimes to ensure functionality with NaTs
    onsets = pd.to_datetime(onset_list)
    
    # Construct a mask, which is True on NaT values
    natmask = pd.isnull(onsets)
    
    # Also keep track of the amount of NaTs
    nat_amount = len(np.where(natmask)[0])
    
    # Use the NOT mask to get all the values which are NOT NaTs
    masked_onsets = onsets[~natmask]
    
    # Sort the times in chronological order
    masked_onsets = np.sort(masked_onsets)
    
    # Returns unique values from the masked array, and an array corresponding to the 
    # mount of each unique instance
    uniques, counts = np.unique(masked_onsets, return_counts=True)
    
    # ...Again convert to pandas datetimes
    uniques = pd.to_datetime(uniques)
    
    # Add the amount of NaTs to the end of counts
    counts = np.append(counts, nat_amount)

    count_sum = sum(counts)
    
    if prints:
        for i, uni_val in enumerate(uniques):
            print(f"{str(uni_val)[:-4]},  n={counts[i]}  ({np.round(100*counts[i]/count_sum, 2)} %)")

        print(f"NaTs: {counts[-1]}  ({np.round(100*counts[-1]/count_sum, 2)} %)")
        print("#----------------------------------")
    
    # Check that not ALL onsets are NaTs (this is evident if counts only holds one number)
    if len(counts) == 1:
        print(f" function unique_vals() output: {onsets[0]}, {counts}")
        return [onsets[0]], counts
    else:
        return uniques, counts

#================================================================================================================

def calc_chnl_nominal_energy(e_min, e_max, mode='gmean'):
    '''
    Calculates the nominal energies of each channel, given the lower and upper bound of
    the channel energies. Channel energy should be given in eVs.
    
    Returns:
        energies: np.array containing the nominal channel energy for each channel
    '''

    #geometric mean
    if mode=='gmean':
        energies = np.sqrt(np.multiply(e_min,e_max))
        return energies
    
    #upper bound of energy channels
    elif mode=='upper':
        energies = e_max
        return energies
    
    else:
        raise Exception("Choose either 'gmean' or 'upper' as the nominal energy.")

#==============================================================================

def sc_distance(mag_data):

    from sunpy.coordinates import frames, get_horizons_coord

    pos = get_horizons_coord("Solar Orbiter", mag_data.index, "id")  # (lon, lat, radius) in (deg, deg, AU)
    pos = pos.transform_to(frames.HeliographicCarrington)
    dist = np.round(pos.radius.value, 2)

    return dist

#==============================================================================

def seek_fit_and_errors(x,y,xerr,yerr, guess=None):
    '''
    Seeks a fit and its errors by ODR (Orthogonal Distance Regression) 
    https://www.tutorialspoint.com/scipy/scipy_odr.htm

    Parameters
    -----------
    x, y :  arrays
            values for the x and y axes
    xerr, yerr : array
            errors of x and y, respectively
    guess : {array-like, list, tuple} with len()==2, default None
            The first guess values for the ODR fit.

    Returns:
    ----------
    out : {ODR.regression}
            An object containing information of the fit and it's errors
    '''

    from scipy.odr import ODR, Model, RealData

    # Define a polynomial of rank 1 to fit the data with.
    def linear_func(p, x):
        slope, const = p
        return slope*x + const

    # Create a model for fitting.
    linear_model = Model(linear_func)

    # Create a RealData object using our initiated data from above.
    data = RealData(x, y, xerr, yerr)

    if not guess:
        slope0, const0 = 1.0, y[-1]
    else:
        slope0, const0 = guess[0], datetime_to_sec([pd.to_datetime(guess[1])])[0]

    # Set up ODR with the model and data.
    odr = ODR(data, linear_model, beta0=[slope0, const0])

    # Run the regression.
    out = odr.run()

    # The built-in pprint method can be used to show the results
    # out.pprint()

    return out


def calculate_cusum_window(time_reso, window_minutes:int=30) -> int:
    """
    Calculates the cusum window in terms of datapoints.
    
    Parameters:
    -----------
    time_reso: str or None
                A pandas-compatible time resolution string. Examples: '20s' or '45min'
    window_minutes: int, default 30
                The amount of minutes of data that a cusum window corresponds to.
    
    Returns:
    ----------
    cusum_window: int
                Cusum window in terms of datapoints.
    """

    if isinstance(time_reso, (pd._libs.tslibs.offsets.Second, pd._libs.tslibs.offsets.Minute, pd._libs.tslibs.offsets.Hour)):
        time_reso = time_reso.freqstr

    if time_reso[-3:] == "min":
        datapoint_multiplier = 1
        reso_value = float(time_reso[:-3]) if len(time_reso) > 3 else 1
    elif time_reso[-1] == 'T':
        datapoint_multiplier = 1
        reso_value = int(time_reso[:-1]) if len(time_reso) > 1 else 1
    elif time_reso[-1] == 'H':
        datapoint_multiplier = 1/60
        reso_value = float(time_reso[:-1]) if len(time_reso) > 1 else 1
    elif time_reso[-1] in ('S', 's'):
        datapoint_multiplier = 60
        reso_value = int(time_reso[:-1]) if len(time_reso) > 1 else 1
    else:
        raise Exception(f"Time resolution format ({time_reso}) not recognized. Use either 'min' or 's'.")

    cusum_window = (window_minutes*datapoint_multiplier)/reso_value
    
    return int(cusum_window)



def set_fig_ylimits(ax:plt.Axes, ylim:list=None, flux_series:pd.Series=None):
    """
    Sets the vertical axis limits of the figure, given an Axes and a tuple or list of y-values.
    """

    # In case not otherwise specified, set the lower limit to half of the smallest plotted value,
    # and the higher limit to 1.5 times the highest value
    if ylim is None:
        ylim = [np.nanmin(flux_series[flux_series > 0]) * 0.5,
                np.nanmax(flux_series) * 1.5]

        # If the lower y-boundary is some ridiculously small number, adjust the y-axis a little
        if np.log10(ylim[1])-np.log10(ylim[0]) > 10:
            ylim[0] = 1e-4

    ylim = ax.set_ylim(ylim)

    return ylim


def set_standard_ticks(ax):
    """
    Handles tickmarks, their sizes etc...
    """
    ticklen = 11
    tickw = 2.8

    ax.tick_params(which="major", length=ticklen, width=tickw, labelsize=TICK_LABELSIZE)
    ax.tick_params(which="minor", length=ticklen-3, width=tickw-0.6)


def set_legend(ax: plt.Axes, legend_loc: str, fontsize:int):

    # Legend placement:
    if legend_loc=="out":
        # loc=3 means that the legend handle is "lower left"
        legend_handle, legend_bbox = 3, (1.0, 0.01)
    elif legend_loc=="in":
        # loc=4 means that the legend handle is "lower right"
        legend_handle, legend_bbox = 4, (1.0, 0.01)
    else:
        raise ValueError(f"Argument legend_loc has to be either 'in' or 'out', not {legend_loc}")

    # Sets the legend
    ax.legend(loc=legend_handle, bbox_to_anchor=legend_bbox, fontsize=fontsize)


def get_figdate(dt_array):
    """
    Gets the string representation as YYYY-MM-DD from a list of datetimes, ignoring NaTs.
    """

    figdate = pd.NaT

    for date in dt_array:

        if not pd.isnull(date):
            figdate = date
            break

    try:
        date = figdate.date().strftime("%Y-%m-%d")

    # An attributeerror is cause by figdate being numpy.datetime64 -object. handle
    # it appropriately
    except AttributeError:
        date = np.datetime_as_string(figdate, unit='D')

    return date


def _isnotebook():
    # https://stackoverflow.com/a/39662359/2336056
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False

# ============================================================================================================================
# this will be ran when importing
# set Jupyter notebook cells to 100% screen size:
if _isnotebook():
    from IPython.core.display import HTML, display
    display(HTML(data="""<style> div#notebook-container { width: 99%; } div#menubar-container { width: 85%; } div#maintoolbar-container { width: 99%; } </style>"""))
