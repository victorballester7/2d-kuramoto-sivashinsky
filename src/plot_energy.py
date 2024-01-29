import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import numpy as np
import sys
import time
import warnings
# import read_data.py
from read_data import read_data_energy, read_data_energy_return
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
warnings.filterwarnings("error")  # treat warnings as errors

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


def plot_energy(filename_E: str, filename_E_return: str,
                t_min: float, save: bool) -> None:
    t = []
    E = []
    u_pi_pi = []
    dE = []
    En = []
    En_1 = []
    x_max = []
    y_max = []

    try:
        data_E = read_data_energy(filename_E)
    except UserWarning:
        data_E = [[], [], [], [], [], []]
    else:
        try:
            idx_1 = np.where(data_E[:, 0] > t_min)[0][0]
        except IndexError:
            idx_1 = 0
        try:
            t = data_E[idx_1:, 0]
            E = data_E[idx_1:, 1]
            dE = data_E[idx_1:, 2]
            u_pi_pi = data_E[idx_1:, 3]
            x_max = data_E[idx_1:, 4]
            y_max = data_E[idx_1:, 5]
        except IndexError:
            t = []
            E = []
            dE = []
            u_pi_pi = []
            x_max = []
            y_max = []

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
    fig, ax = plt.subplots(1, 5, figsize=(12, 6))
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

    ax[3].plot(t, u_pi_pi, 'g-')
    ax[3].set_xlabel(r"$t$")
    ax[3].set_ylabel(r'$u(\pi, \pi)$')

    ax[4].plot(x_max, y_max, 'r.')
    ax[4].set_xlabel(r"$x_{max}$")
    ax[4].set_ylabel(r'$y_{max}$')
    # ax[4].set_aspect('equal')
    ax[4].set_xlim(0, 2 * np.pi)
    ax[4].set_ylim(0, 2 * np.pi)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.plot(t, E, '-')
    # data_E2 = read_data_energy("data/energy2.txt")
    # try:
    #     idx_1 = np.where(data_E2[:, 0] > t_min)[0][0]
    # except IndexError:
    #     idx_1 = 0
    # try:
    #     t2 = data_E2[idx_1:, 0]
    #     u2_pi_pi = data_E2[idx_1:, 3]
    # except IndexError:
    #     t2 = []
    #     u2_pi_pi = []
    # increase font size
    # I want to select the x-range for the zoomed region. I have figured it out suitable values
    # by trial and error. How can I pass more elegantly the dates as something
    # like
    # x1 = 205
    # x2 = 310

    # # select y-range for zoomed region
    # y1 = -50
    # y2 = 50
    # ax.set_ylim(-1000, 500)
    # Make the zoom-in plot:
    font_size = 40
    plt.rcParams.update({'font.size': font_size})
    # # increase size xticks and yticks
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlabel(r"$t$", fontsize=font_size)
    plt.ylabel(r'$E$', fontsize=font_size)

    # axins = zoomed_inset_axes(
    #     ax,
    #     6,
    #     loc='lower left',
    #     bbox_to_anchor=(
    #         -0.66,
    #         0.05),
    #     bbox_transform=ax.transAxes)

    # axins.plot(E, dE, '-')
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
    # axins.set_aspect(0.75)
    # mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")
    # plt.xticks(fontsize=font_size // 1.5)
    # plt.yticks(fontsize=font_size // 1.5)

    # save the figure
    plt.savefig("latex/images/tw_energy.pdf", bbox_inches='tight')

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
