import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from matplotlib import colormaps as cm
from matplotlib.ticker import LinearLocator
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes
from matplotlib.text import Text
from typing import cast
from abc import ABC, abstractmethod
from typing import Any
from matplotlib.artist import Artist
from matplotlib.colors import Normalize
import sys
import time
# incl

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


def read_data(filename: str) -> np.ndarray:
    """Read data from filename. The data is of the form:
    t_1 E_1 dE_1
    t_2 E_2 dE_2
    ...

    Args:
        filename (str): File path.

    Returns:
        np.ndarray: Data.
    """
    data = np.loadtxt(filename)
    return data


def plot_energy(filename: str) -> None:
    data = read_data(filename)
    # get the first index of t such that t > 100
    # idx = np.where(data[:, 0] > 100)[0][0]
    idx = 0
    t, E, dE = data[idx:, 0], data[idx:, 1], data[idx:, 2]
    # plot t-E and next to it E-dE
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    ax[0].plot(t, E, '-')
    ax[0].set_xlabel(r"$t$")
    ax[0].set_ylabel(r'$E$')
    ax[1].plot(E, dE, '-')
    ax[1].set_xlabel(r"$E$")
    # label dE as \dot{E}
    ax[1].set_ylabel(r'$\dot{E}$')
    # ax[1].set_aspect('equal')
    # add space between subplots
    fig.subplots_adjust(wspace=0.3)
    end_time = time.time()
    print("Total time for plotting: ", int(
        UNIT_TIME*(end_time - start_time)), LABEL_TIME)
    plt.show()


UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
filename = "data/energy.txt"
plot_energy(filename)
