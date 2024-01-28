import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import numpy as np
import sys
import time
import warnings
from matplotlib import colormaps as cm
from typing import cast
from mpl_toolkits.mplot3d import Axes3D
from read_data import read_data_solution
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings("error")  # treat warnings as errors

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


def scale_text(fig, scale_factor):
    """
    Scale the text in Matplotlib figure by a scale factor.

    Parameters:
        - fig: Matplotlib figure object
        - scale_factor: Scale factor for the text (e.g., 1.2 for 20% larger)
    """
    # Iterate through all the text elements in the figure
    for text in fig.findobj(plt.Text):
        # Get the current font size
        current_font_size = text.get_fontsize()
        # Scale the font size
        new_font_size = current_font_size * scale_factor
        # Set the new font size
        text.set_fontsize(new_font_size)


def plot_solution(filename: str, nu1: float, nu2: float,
                  time_plot: float) -> None:
    Z = read_data_solution(filename)
    # Z is of the form
    # k_00 k_01 ... k_0n
    # k_10 k_11 ... k_1n
    # ...
    # k_n0 k_n1 ... k_nn
    MAX_FREQ = 16
    Z = Z[:MAX_FREQ, :MAX_FREQ]
    Z = Z / np.max(Z)
    Z_MIN = np.min(Z)
    Z_MAX = np.max(Z)
    # create a colormap from white to magenta
    colors = [(1, 1, 1), (1, 0, 0), (0.6, 0, 0), (0.2, 0, 0)]  # White to Red

    # Create the custom colormap
    custom_cmap = LinearSegmentedColormap.from_list(
        "Custom White to Red", colors)

    color = 'ocean'
    X_LABEL = r'$k_x$'
    Y_LABEL = r'$k_y$'

    font_size = 40

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(Z, cmap=custom_cmap, vmin=Z_MIN, vmax=Z_MAX)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_xticks(np.arange(0, MAX_FREQ, 4))
    ax.set_yticks(np.arange(0, MAX_FREQ, 4))
    ax.invert_yaxis()  # invert y axis
    # add colorbar to the right
    cbar = ax.figure.colorbar(ax.images[0], ax=ax, orientation='vertical')

    # Scale the text in the current figure by a scale factor of 1.2 (20%
    # larger)
    scale_text(plt.gcf(), 1.5)
    # save figure
    bbox = fig.bbox_inches.from_bounds(0.49, -0.01, 5.5, 4.5)
    fig.savefig(
        "latex/images/slice_freq_nu1_" +
        str(nu1) +
        "_nu2_" +
        str(nu2) +
        "_time_" +
        str(time_plot) +
        ".pdf", bbox_inches=bbox)
    end_time = time.time()
    # plt.show()
    print("Total time for plotting: ", int(
        UNIT_TIME * (end_time - start_time)), LABEL_TIME)


UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
filename = "data/freq_slice.txt"
try:
    nu1 = float(sys.argv[1])
    nu2 = float(sys.argv[2])
    time_plot = float(sys.argv[3])
except IndexError:
    nu1 = np.nan
    nu2 = np.nan
    time_plot = np.nan
plot_solution(filename, nu1, nu2, time_plot)
