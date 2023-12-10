import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
from typing import Union, Any
from matplotlib.cm import ScalarMappable
from matplotlib.artist import Artist
import sys


class Plot(ABC):
    def __init__(self):
        pass

    def read_data_file(self, file_path):
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

    def plot_setup(self):
        # Create a figure
        self.fig = plt.figure()
        self.ax = self.get_ax()

        # Read data
        script_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = script_dir + '/../data/data.txt'  # Replace with the actual file path
        self.headers, self.data_blocks = self.read_data_file(file_path)

        # Create data
        self.nx = self.data_blocks.shape[1]
        self.ny = self.data_blocks.shape[2]

        self.X = np.linspace(0, 2*np.pi, self.nx)
        self.Y = np.linspace(0, 2*np.pi, self.ny)
        self.X, self.Y = np.meshgrid(self.X, self.Y, indexing='ij')
        self.Z = self.data_blocks[0]

        # DEFAULTS
        self.FONTSIZE_TITLE = 15
        self.FONTSIZE_TIME = 10
        self.color = cm['viridis']
        self.FPS = 50
        self.Z_MAX = np.max(np.abs(self.data_blocks))
        self.Z_MIN = -self.Z_MAX
        self.plot_args = self.get_plot_args()
        self.y_pos_text = 1.03
        self.y_pos_title = self.y_pos_text + 0.04
        if self.is_3d():
            self.y_pos_text = 0.95
        self.text_args = {'x': 0.5, 'y': self.y_pos_text, 's': '', 'transform': self.ax.transAxes,
                          'fontsize': self.FONTSIZE_TIME, 'horizontalalignment': 'center'}
        self.time_text = self.get_text()

        # Title
        self.ax.set_title('2D Kuramoto-Sivashinsky equation',
                          fontsize=self.FONTSIZE_TITLE, fontweight='bold', y=self.y_pos_title)
        # Axes
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        if self.is_3d():
            self.ax = cast(Axes3D, self.ax)
            self.ax.set_zlabel('z')
            self.ax.set_zlim(self.Z_MIN, self.Z_MAX)

        # Create plot
        self.plot = [self.get_plot()]

        # Add a color bar which maps values to colors.
        self.fig.colorbar(self.plot[0], shrink=0.5, aspect=10, location='left')

        # Animation update function
        def init():
            return self.plot[0], self.time_text

        def update(frame):
            # Clear the previous frame
            # ax.clear()
            self.plot[0].remove()

            # Update the arrays
            self.Z = self.data_blocks[frame]
            # Z = np.sin(2*X) + np.cos(2*Y+10*headers[frame])

            # Update the time text
            self.time_text.set_text('t = %.2f' % self.headers[frame])

            # Update the plot
            self.plot[0] = self.get_plot()

            return self.plot[0], self.time_text

        # Create the animation
        ani = FuncAnimation(self.fig, update, frames=len(self.data_blocks),
                            interval=1000/self.FPS, blit=False, init_func=init)
        # Show the animation
        plt.show()

    def get_ax(self) -> Axes | Axes3D:
        if self.is_3d():
            return cast(Axes3D, self.fig.add_subplot(projection='3d'))
        else:
            return self.fig.add_subplot()

    def get_text(self) -> Text:
        if self.is_3d():
            self.ax = cast(Axes3D, self.ax)
            return self.ax.text2D(**self.text_args)
        else:
            return self.ax.text(**self.text_args)

    @abstractmethod
    def is_3d(self) -> bool:
        pass

    @abstractmethod
    def get_plot(self) -> Artist:
        pass

    @abstractmethod
    def get_plot_args(self) -> dict[str, Any]:
        pass


class Surface(Plot):
    def __init__(self):
        super().__init__()
        self.plot_setup()

    def is_3d(self) -> bool:
        return True

    def get_plot(self):
        # assume self.ax is Axes3D
        self.ax = cast(Axes3D, self.ax)
        return self.ax.plot_surface(self.X, self.Y, self.Z, **self.plot_args)

    def get_plot_args(self) -> dict[str, Any]:
        return {'rstride': 1, 'cstride': 1, 'cmap': self.color, 'linewidth': 0.01,
                'antialiased': True, 'color': 'w', 'shade': True, 'vmin': self.Z_MIN, 'vmax': self.Z_MAX}


class Contour(Plot):
    def __init__(self):
        super().__init__()
        self.plot_setup()

    def is_3d(self) -> bool:
        return False

    def get_plot(self):
        return self.ax.contourf(self.X, self.Y, self.Z, **self.plot_args)

    def get_plot_args(self) -> dict[str, Any]:
        plt.style.use('_mpl-gallery-nogrid')
        self.num_levels = 20
        self.levels = np.linspace(self.Z_MIN, self.Z_MAX, self.num_levels)
        return {'cmap': self.color, 'levels': self.levels}


arg = sys.argv[1]
if int(arg) == 1:
    Surface()
else:
    Contour()
# plot_surface()
# plot_contour()
