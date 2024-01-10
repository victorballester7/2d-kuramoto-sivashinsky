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


def read_data_energy(filename: str) -> np.ndarray:
    """Read data from filename. The data is of the form:
    t_1 E_1 dE_1 u(pi, pi)
    t_2 E_2 dE_2 u(pi, pi)
    ...

    Args:
        filename (str): File path.

    Returns:
        np.ndarray: Data.
    """
    data = np.loadtxt(filename)
    return data


def read_data_energy_return(filename: str) -> np.ndarray:
    """Read data from filename. The data is of the form:
    t1 E_1
    t2 E_2
    ...

    Args:
        filename (str): File path.

    Returns:
        np.ndarray: Data.
    """
    data = np.loadtxt(filename)
    return data


def plot_energy(filename_E: str, filename_E_return: str, t_min: float) -> None:
    data_E = read_data_energy(filename_E)
    data_E_return = read_data_energy_return(filename_E_return)

    # get the first index of t such that t > 100
    try:
        idx_1 = np.where(data_E[:, 0] > t_min)[0][0]
        idx_2 = np.where(data_E_return[:, 0] > t_min)[0][0]
    except IndexError:
        idx_1 = 0
        idx_2 = 0
    En = data_E_return[idx_2:-1, 1]
    En_1 = data_E_return[idx_2 + 1:, 1]

    # idx = 0
    t, E, dE = data_E[idx_1:, 0], data_E[idx_1:, 1], data_E[idx_1:, 2]
    # plot t-E and next to it E-dE
    fig, ax = plt.subplots(1, 3, figsize=(12, 6))
    ax[0].plot(t, E, '-')
    ax[0].set_xlabel(r"$t$")
    ax[0].set_ylabel(r'$E$')

    ax[1].plot(E, dE, '-')
    ax[1].set_xlabel(r"$E$")
    # label dE as \dot{E}
    ax[1].set_ylabel(r'$\dot{E}$')
    # ax[1].set_aspect('equal')

    ax[2].plot(En, En_1, '.')
    ax[2].set_xlabel(r"$E_n$")
    ax[2].set_ylabel(r'$E_{n+1}$')

    # add space between subplots
    fig.subplots_adjust(wspace=0.5)
    end_time = time.time()
    print("Total time for plotting: ", int(
        UNIT_TIME*(end_time - start_time)), LABEL_TIME)
    plt.show()


UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
filename_E = "data/energy.txt"
filename_E_return = "data/energy_return.txt"
try:
    t_min = float(sys.argv[1])
except IndexError:
    t_min = 200.0
plot_energy(filename_E, filename_E_return, t_min)
