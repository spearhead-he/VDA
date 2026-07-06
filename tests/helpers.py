# tests/helpers.py
def strip_figure_text(fig):
    """
    Remove all text elements that vary across matplotlib versions.
    Use before returning a figure in mpl_image_compare tests.
    """
    fig.suptitle('')
    for text in fig.texts:
        text.set_text('')
    for ax in fig.axes:
        ax.set_title('')
        ax.set_xlabel('')      # not removed by remove_text=True
        ax.set_ylabel('')      # not removed by remove_text=True
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        leg = ax.get_legend()
        if leg:
            leg.remove()
    return fig