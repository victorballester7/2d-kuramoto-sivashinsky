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
    Z_MIN = np.min(Z)
    Z_MAX = np.max(Z)
    # Z_MAX = 7
    # Z_MIN = -7
    color = cm['inferno']
    X_LABEL = 'x'
    Y_LABEL = 'y'
    Z_LABEL = 'u'

    nx = Z.shape[0]
    ny = Z.shape[1]
    font_size = 40
    X = np.linspace(0, 2 * np.pi, nx)
    Y = np.linspace(0, 2 * np.pi, ny)

    X, Y = np.meshgrid(X, Y, indexing='ij')

    fig = plt.figure()
    ax = cast(Axes3D, fig.add_subplot(111, projection='3d'))
    # change color lines to light gray
    plot_args = {'rstride': 1, 'cstride': 1, 'cmap': color, 'linewidth': 0.01,
                 'antialiased': True, 'shade': True, 'vmin': Z_MIN, 'vmax': Z_MAX}
    num_levels = 14
    offset_rate = 2
    levels = np.linspace(Z_MIN, Z_MAX, num_levels)
    plot_args_extra = {'zdir': 'z', 'offset': offset_rate *
                       Z_MIN, 'cmap': color, 'levels': levels}
    ax.set_zlim(offset_rate * Z_MIN, Z_MAX)
    ax.plot_surface(X, Y, Z, **plot_args)
    ax.contour(X, Y, Z, **plot_args_extra)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_zlabel(Z_LABEL)
    # Scale the text in the current figure by a scale factor of 1.2 (20%
    # larger)
    scale_text(plt.gcf(), 1.5)
    # save figure
    bbox = fig.bbox_inches.from_bounds(1.55, 0.44, 4, 3.6)
    fig.savefig(
        "latex/images/slice_nu1_" +
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
filename = "data/solution_slice.txt"
try:
    nu1 = float(sys.argv[1])
    nu2 = float(sys.argv[2])
    time_plot = float(sys.argv[3])
except IndexError:
    nu1 = np.nan
    nu2 = np.nan
    time_plot = np.nan
plot_solution(filename, nu1, nu2, time_plot)
