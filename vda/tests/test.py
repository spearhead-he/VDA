import filecmp
import math
import matplotlib
import pytest

from warnings import simplefilter, filterwarnings
from pandas.errors import PerformanceWarning
from astropy.visualization import quantity_support

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
pytest --mpl-generate-path=vda/tests/baseline vda/tests/test.py

To run the tests locally, go to the base directory of the repository and run:
pytest -rP --mpl --mpl-baseline-path=baseline --mpl-baseline-relative --mpl-generate-summary=html --cov=vda vda/tests/test.py
"""


@pytest.mark.mpl_image_compare(remove_text=True, deterministic=True)
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

    # particle selection (not needed)
    vda_displayer.display_particle_selection()

    # load data
    vda.construct_particles_df()

    vda.construct_energies_df()

    vda_displayer.display_groupings()

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

    return fig
