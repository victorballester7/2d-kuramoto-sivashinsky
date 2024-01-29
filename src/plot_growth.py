import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import numpy as np
import sys
import time
import warnings
# import read_data.py
from read_data import read_data_energy, read_data_energy_return
warnings.filterwarnings("error")  # treat warnings as errors

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


def plot_growth(filename: str) -> None:
    t = []
    u_max = []

    u_max = read_data_energy_return(filename)

    t = u_max[:, 0]
    u_max = u_max[:, 1]
    # Now you can use t, E, dE, En, and En_1 in the rest of your code

    # plot t-E and next to it E-dE
    font_size = 40
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.plot(t, u_max, 'g-')
    ax.set_xlabel(r"$t$", fontsize=font_size)
    ax.set_ylabel('Relative error', fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)

    # fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    # ax.plot(t, E, '-')
    # # data_E2 = read_data_energy("data/energy2.txt")
    # # try:
    # #     idx_1 = np.where(data_E2[:, 0] > t_min)[0][0]
    # # except IndexError:
    # #     idx_1 = 0
    # # try:
    # #     t2 = data_E2[idx_1:, 0]
    # #     u2_pi_pi = data_E2[idx_1:, 3]
    # # except IndexError:
    # #     t2 = []
    # #     u2_pi_pi = []
    # # increase font size
    # # I want to select the x-range for the zoomed region. I have figured it out suitable values
    # # by trial and error. How can I pass more elegantly the dates as something
    # # like
    # # x1 = 205
    # # x2 = 310

    # # # select y-range for zoomed region
    # # y1 = -50
    # # y2 = 50
    # # ax.set_ylim(-1000, 500)
    # # Make the zoom-in plot:
    # font_size = 40
    # plt.rcParams.update({'font.size': font_size})
    # # # increase size xticks and yticks
    # plt.xticks(fontsize=font_size)
    # plt.yticks(fontsize=font_size)
    # plt.xlabel(r"$t$", fontsize=font_size)
    # plt.ylabel(r'$E$', fontsize=font_size)

    # # axins = zoomed_inset_axes(
    # #     ax,
    # #     6,
    # #     loc='lower left',
    # #     bbox_to_anchor=(
    # #         -0.66,
    # #         0.05),
    # #     bbox_transform=ax.transAxes)

    # # axins.plot(E, dE, '-')
    # # axins.set_xlim(x1, x2)
    # # axins.set_ylim(y1, y2)
    # # axins.set_aspect(0.75)
    # # mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")
    # # plt.xticks(fontsize=font_size // 1.5)
    # # plt.yticks(fontsize=font_size // 1.5)

    # # save the figure
    # plt.savefig("latex/images/hetero_chaotic.pdf", bbox_inches='tight')

    # add space between subplots
    end_time = time.time()
    print("Total time for plotting: ", int(
        UNIT_TIME * (end_time - start_time)), LABEL_TIME)
    # write nu1, nu2 with 2 decimal places
    plt.savefig("latex/images/burst_nu1={}_nu2={}.pdf".format(
        "{:.2f}".format(nu1), "{:.2f}".format(nu2)), bbox_inches='tight')


UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
# filename_E = "others/2D-Kuramoto-Sivashinsky/u_norm_nu1=0-15_nu2=0-03.txt"
filename_growth = "data/growth_rate_u.txt"

try:
    nu1 = float(sys.argv[1])
    nu2 = float(sys.argv[2])
except IndexError:
    nu1 = np.nan
    nu2 = np.nan
plot_growth(filename_growth)
