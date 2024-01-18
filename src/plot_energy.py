import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import numpy as np
import sys
import time
import warnings

warnings.filterwarnings("error")  # treat warnings as errors

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


def plot_energy(filename_E: str, filename_E_return: str,
                t_min: float, save: bool) -> None:
    t = []
    E = []
    u_pi_pi = []
    dE = []
    En = []
    En_1 = []

    try:
        data_E = read_data_energy(filename_E)
    except UserWarning:
        data_E = [[], [], []]
    else:
        try:
            idx_1 = np.where(data_E[:, 0] > t_min)[0][0]
        except IndexError:
            idx_1 = 0
        try:
            t, E, dE, u_pi_pi = data_E[idx_1:, 0], data_E[idx_1:,
                                                          1], data_E[idx_1:, 2], data_E[idx_1:, 3]
        except IndexError:
            t = []
            E = []
            dE = []
            u_pi_pi = []

    try:
        data_E_return = read_data_energy_return(filename_E_return)
    except UserWarning:
        data_E_return = [[], []]
    else:
        try:
            idx_2 = np.where(data_E_return[:, 0] > t_min)[0][0]
        except IndexError:
            idx_2 = 0
        try:
            En = data_E_return[idx_2:-1, 1]
            En_1 = data_E_return[idx_2 + 1:, 1]
        except IndexError:
            En = []
            En_1 = []

    # Now you can use t, E, dE, En, and En_1 in the rest of your code

    # plot t-E and next to it E-dE
    fig, ax = plt.subplots(1, 4, figsize=(12, 6))
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

    ax[3].plot(t, u_pi_pi, '-')
    ax[3].set_xlabel(r"$t$")
    ax[3].set_ylabel(r'$u(\pi, \pi)$')

    # add space between subplots
    fig.subplots_adjust(wspace=0.5)
    end_time = time.time()
    print("Total time for plotting: ", int(
        UNIT_TIME * (end_time - start_time)), LABEL_TIME)
    if save:
        title = r"$\nu_1 = {}, \nu_2 = {}$".format(nu1, nu2)
        # add the title above the second subplot
        fig.suptitle(title, fontsize=16)
        # write nu1, nu2 with 2 decimal places
        plt.savefig("img/energy/energy_nu1={}_nu2={}.jpg".format(
            "{:.2f}".format(nu1), "{:.2f}".format(nu2)))
    else:
        plt.show()


UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
filename_E = "data/energy.txt"
# filename_E = "others/2D-Kuramoto-Sivashinsky/u_norm_nu1=0-15_nu2=0-03.txt"
filename_E_return = "data/energy_return.txt"
try:
    t_min = float(sys.argv[1])
except IndexError:
    t_min = 200.0

try:
    nu1 = float(sys.argv[2]) / 100.0
    nu2 = float(sys.argv[3]) / 100.0
except IndexError:
    plot_energy(filename_E, filename_E_return, t_min, False)
else:
    plot_energy(filename_E, filename_E_return, t_min, True)
