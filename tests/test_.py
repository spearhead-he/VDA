import filecmp
import math
import matplotlib
import pytest

from warnings import simplefilter, filterwarnings
from pandas.errors import PerformanceWarning
from astropy.visualization import quantity_support

from tests.helpers import strip_figure_text
from vda_tool_configuration import VDA_parameters
from vda_views import VDA_nb_displayer
from vda import VDA

# omit Pandas' PerformanceWarning
simplefilter(action='ignore', category=PerformanceWarning)
filterwarnings(action='ignore', message="Discarding nonzero nanoseconds in conversion")


"""
Install dependencies for tests:
pip install flake8 pytest pytest-doctestplus pytest-cov pytest-mpl

To create/update the baseline images, run the following command from the base package dir:
pytest --mpl-generate-path=tests/baseline tests/test_.py

To run the tests locally, go to the base directory of the repository and run:
pytest -rP --mpl --mpl-baseline-path=baseline --mpl-baseline-relative --mpl-generate-summary=html tests/test_.py
"""

# skip image comparison tests for matplotlib < 3.11
_mpl_old = (int(matplotlib.__version__.split(".")[1])) < 11
_image_compare = pytest.mark.mpl_image_compare(remove_text=True, deterministic=True) if not _mpl_old else lambda f: f


@_image_compare
def test_vda_default():
    vda_parameters = VDA_parameters()
    vda = VDA(vda_parameters)
    vda_displayer = VDA_nb_displayer(vda)

    # not needed
    vda_displayer.display_view_toggle()
    vda_displayer.display_date_range()

    # construct times
    vda.construct_times_df()

    # not needed
    vda_displayer.display_load_data_option()

    # load data
    vda_displayer.construct_energies_df()

    # particle selection (not needed)
    vda_displayer.display_particle_selection()

    vda.construct_particles_df()

    vda.group_energy_channels()

    vda_displayer.display_onset_method_selection()

    vda_displayer.display_onset_method_parameters()

    vda.plot_bg_selection()

    vda.calculate_onsets()

    vda.clean_onsets()

    vda.construct_options_df()

    vda_displayer.display_onset_selection_selection()

    vda_displayer.select_onsets()

    vda.parameters.selected_onsets

    vda.construct_energy_channels_characteristics()

    quantity_support()

    vda.define_spacecraft_parameters()

    fig = vda.plot(savefig=False, returnfig=True)

    # check legend contents manually
    handles, labels = fig.axes[0].get_legend_handles_labels()
    assert labels == ['Linear Regression',
                      'Extra Time = 0:06:40',
                      'Release Time = 2021-10-28 15:31:11 +/- 0:03:38',
                      'APL = 1.76 +/- 0.12']

    if not _mpl_old:
        # Strip before returning — don't rely solely on remove_text=True
        return strip_figure_text(fig)
