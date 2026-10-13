"""End-to-end test of the default event (network test): downloads the Solar Orbiter data and compares the VDA plot with
its baseline image. The commands to run it and to update the baseline image are in docs/development.md."""
import matplotlib
import pytest

from warnings import simplefilter, filterwarnings
from pandas.errors import PerformanceWarning
from astropy.visualization import quantity_support

from tests.helpers import strip_figure_text
from spearhead.vda import VDA, VDA_parameters, views

# omit Pandas' PerformanceWarning
simplefilter(action='ignore', category=PerformanceWarning)
filterwarnings(action='ignore', message="Discarding nonzero nanoseconds in conversion")


# skip image comparison tests for matplotlib < 3.11
_mpl_old = (int(matplotlib.__version__.split(".")[1])) < 11
_image_compare = pytest.mark.mpl_image_compare(remove_text=True, deterministic=True) if not _mpl_old else lambda f: f


@pytest.mark.network
@_image_compare
def test_vda_default():
    vda_parameters = VDA_parameters()
    vda = VDA(vda_parameters)

    vda.construct_times_df()

    vda.construct_energies_df()

    vda.construct_particles_df()

    vda.group_energy_channels()

    for event_no in vda.df_grouped.index.unique(level=0):
        views.plot_bg(vda, event_no)

    vda.calculate_onsets()

    vda.clean_onsets()

    vda.construct_options_df()

    vda.select_onsets()

    vda.construct_energy_channels_characteristics()

    quantity_support()

    vda.define_spacecraft_parameters()

    vda.compute_vda()

    fig = views.plot_vda(vda, 1)

    # check legend contents manually
    handles, labels = fig.axes[0].get_legend_handles_labels()
    assert labels == ['Protons',
                      'Electrons',
                      'Linear Regression',
                      'Extra Time = 0:06:40',
                      'Release Time = 2021-10-28 15:26:11 +/- 0:03:38',
                      'APL = 1.76 +/- 0.12']

    if not _mpl_old:
        # Strip before returning — don't rely solely on remove_text=True
        return strip_figure_text(fig)
