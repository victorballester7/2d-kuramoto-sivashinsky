import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import os


def read_data_file(file_path):
    headers = []
    data_blocks = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line:
            header = float(line)
            headers.append(header)

            i += 1
            block_lines = []
            while i < len(lines) and lines[i].strip():
                aux = list(map(float, lines[i].split()))
                # repeat the first element of each line to make the animation periodic
                aux.append(aux[0])
                block_lines.append(aux)
                i += 1
            # repeat the first line to make the animation periodic
            block_lines.append(block_lines[0])
            data_block = np.array(block_lines)
            data_blocks.append(data_block)
        i += 1

    headers_array = np.array(headers)
    data_blocks_array = np.array(data_blocks)

    # now we dupplicate the

    return headers_array, data_blocks_array


def plot():
    # Create a figure and a 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # cast ax to Axes3D
    # ax = Axes3D(fig)

    # Read data
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = script_dir + '/../data/data.txt'  # Replace with the actual file path
    # Replace with the actual file path
    file_path_sol = script_dir + '/../data/sol_aux.txt'
    headers, data_blocks = read_data_file(file_path)
    # headers_sol, data_blocks_sol = read_data_file(file_path_sol)
    # for i in range(len(data_blocks)):
    #     print(np.max(np.abs(data_blocks[i] - data_blocks_sol[i])))

    # Make data.
    nx = data_blocks.shape[1]
    ny = data_blocks.shape[2]

    X = np.linspace(0, 2*np.pi, nx)
    Y = np.linspace(0, 2*np.pi, ny)
    X, Y = np.meshgrid(X, Y, indexing='ij')
    DATA_PLOT = data_blocks
    Z = DATA_PLOT[0]
    # Z = np.sin(2*Y) + np.cos(2*Y)

    # DEFAULTS
    FONTSIZE_TITLE = 15
    FONTSIZE_TIME = 10

    # Customize axis
    ax.set_title('2D Kuramoto-Sivashinsky equation',
                 fontsize=FONTSIZE_TITLE, fontweight='bold')
    # center text
    time_text = ax.text2D(0.5, 0.95, '', transform=ax.transAxes,
                          fontsize=FONTSIZE_TIME, horizontalalignment='center')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    Z_MIN = -2
    Z_MAX = 2
    ax.set_zlim(Z_MIN, Z_MAX)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    # ax.zaxis.set_major_formatter('{x:.02f}')

    # Plot the surface.
    plot_args = {'rstride': 1, 'cstride': 1, 'cmap': cm.coolwarm,
                 'linewidth': 0.01, 'antialiased': True, 'color': 'w', 'shade': True, 'vmin': Z_MIN, 'vmax': Z_MAX}
    # surf = [ax.plot_surface(X, Y, Z, **plot_args)]
    surf = ax.plot_surface(X, Y, Z, **plot_args)
    # Add a color bar which maps values to colors.
    cb = fig.colorbar(surf, shrink=0.5, aspect=10,
                      location='left')
    cax = cb.ax

    # Animation update function

    def init():
        return surf, time_text

    def update(frame, data_blocks, surf):
        # Clear the previous frame
        ax.clear()
        # surf.remove()
        # surf[0].remove()
        # Update the arrays
        Z = data_blocks[frame]

        # Update the time text
        time_text.set_text('t = %.2f' % headers[frame])

        # Update the scatter plot
        surf = ax.plot_surface(X, Y, Z, **plot_args)

        # update the colorbar
        # fig.colorbar(surf[0], cax=cax)
        # cb.update_normal(surf[0])

        return surf, time_text

    # Create the animation
    ani = FuncAnimation(fig, update, frames=len(DATA_PLOT), fargs=(
        DATA_PLOT, surf), interval=50, save_count=100, blit=True, init_func=init)
    # ani = FuncAnimation(fig, update, frames=len(
    # data_blocks), fargs=(data_blocks, surf), interval=50, save_count=100, blit=False, init_func=init)

    # Show the animation
    plt.show()


plot()
