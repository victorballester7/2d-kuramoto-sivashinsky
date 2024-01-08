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


class Plot(ABC):
    # DEFAULTS
    TITLE = '2D Kuramoto-Sivashinsky equation'
    X_LABEL = 'x'
    Y_LABEL = 'y'
    Z_LABEL = 'z'
    FONTSIZE_TITLE = 15
    FONTSIZE_TIME = 10
    y_pos_text = 1.03
    y_pos_title = y_pos_text + 0.14
    color = cm['viridis']
    color_extra = 'royalblue'
    colorbar_args = {'shrink': 0.5, 'aspect': 10, 'location': 'left'}
    FPS = 50
    speed = 1

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
        # Replace with the actual file path
        file_path = script_dir + '/../data/solution.txt'
        self.headers, self.data_blocks = self.read_data_file(file_path)

        # Create data
        self.nx = self.data_blocks.shape[1]
        self.ny = self.data_blocks.shape[2]

        self.X = np.linspace(0, 2*np.pi, self.nx)
        self.Y = np.linspace(0, 2*np.pi, self.ny)
        self.X, self.Y = np.meshgrid(self.X, self.Y, indexing='ij')
        self.Z = self.data_blocks[0]

        self.Z_MAX = np.max(np.abs(self.data_blocks))
        self.Z_MIN = -self.Z_MAX
        self.plot_args = self.get_plot_args()
        self.plot_args_extra = self.get_plot_args_extra()
        if self.is_3d():
            self.y_pos_text = 0.95
        self.text_args = {'x': 0.5, 'y': self.y_pos_text, 's': '', 'transform': self.ax.transAxes,
                          'fontsize': self.FONTSIZE_TIME, 'horizontalalignment': 'center'}
        self.time_text = self.get_text()

        # Title
        # plt.title(self.TITLE,
        # fontsize=self.FONTSIZE_TITLE, fontweight='bold', y=self.y_pos_title)

        # subtitle
        # self.subtitle = r'$\nu_1 = ' + \
        #     str(self.nu1) + r'$, $\nu_2 = ' + str(self.nu2) + r'$'
        # self.subtitle_args = {'x': 0.5, 'y': self.y_pos_title - 0.14, 's': self.subtitle, 'transform': self.ax.transAxes,
        #                       'fontsize': self.FONTSIZE_TITLE - 4, 'horizontalalignment': 'center'}
        # plt.title(
        #     subtitle, fontsize=self.FONTSIZE_TITLE - 2, y=self.y_pos_title - 0.05)

        # Axes
        self.ax.set_xlabel(self.X_LABEL)
        self.ax.set_ylabel(self.Y_LABEL)
        if self.is_3d():
            self.ax = cast(Axes3D, self.ax)
            self.ax.set_zlabel(self.Z_LABEL)
            self.ax.set_zlim(self.Z_MIN, self.Z_MAX)
            # self.ax.text2D(**self.subtitle_args)
        # else:
            # self.ax.text(**self.subtitle_args)

        # Create plot
        self.plot = [i for i in self.get_plot()]

        self.num_plots = len(self.plot)

        # Add a color bar which maps values to colors.
        # mappable = self.plot[0] if isinstance(self.plot[0], Artist) else self.plot[0][0]
        # self.fig.colorbar(mappable=m, shrink=0.5, aspect=10, location='left')

        # self.fig.colorbar(self.plot[0], **self.colorbar_args)

        # Animation update function
        def init():
            return self.plot[0], self.time_text

        def update(frame):
            # Clear the previous frame
            # ax.clear()
            for plot in self.plot:
                plot.remove()

            # Update the arrays
            self.Z = self.data_blocks[frame * self.speed]
            # Z = np.sin(2*X) + np.cos(2*Y+10*headers[frame])

            # Update the time text
            self.time_text.set_text('t = %.2f' %
                                    self.headers[frame * self.speed])

            # Update the plot
            # self.plot[0] = self.get_plot()
            tmp = self.get_plot()
            for i in range(self.num_plots):
                self.plot[i] = tmp[i]

            return self.plot[0], self.time_text

        # Create the animation
        self.ani = FuncAnimation(self.fig, update, frames=len(self.data_blocks) // self.speed,
                                 interval=1000/self.FPS, blit=False, init_func=init)

        end_time = time.time()
        print("Total time for animating: ", int(
            UNIT_TIME*(end_time - start_time)), LABEL_TIME)

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

    def show_plot(self):
        plt.show()

    @abstractmethod
    def is_3d(self) -> bool:
        pass

    @abstractmethod
    def get_plot(self) -> list[Any]:
        pass

    @abstractmethod
    def get_plot_args(self) -> dict[str, Any]:
        pass

    @abstractmethod
    def get_plot_args_extra(self) -> dict[str, Any]:
        pass


class Surface(Plot):
    def __init__(self):
        super().__init__()
        self.plot_setup()
        self.show_plot()

    def is_3d(self) -> bool:
        return True

    def get_plot(self):
        # assume self.ax is Axes3D
        self.ax = cast(Axes3D, self.ax)
        return [self.ax.plot_surface(self.X, self.Y, self.Z, **self.plot_args)]

    def get_plot_args(self) -> dict[str, Any]:
        return {'rstride': 1, 'cstride': 1, 'cmap': self.color, 'linewidth': 0.01,
                'antialiased': True, 'color': 'w', 'shade': True, 'vmin': self.Z_MIN, 'vmax': self.Z_MAX}

    def get_plot_args_extra(self):
        pass


class Contour(Plot):
    def __init__(self):
        super().__init__()
        self.plot_setup()
        # self.custom_colorbar()
        self.show_plot()

    # def custom_colorbar(self):
    #     # remove the previous colorbar
    #     self.fig.delaxes(self.fig.axes[1])

    #     # create a new colorbar with self.levels
    #     norm = Normalize(vmin=self.Z_MIN, vmax=self.Z_MAX)
    #     # a previous version of this used
    #     # norm= matplotlib.colors.Normalize(vmin=cs.vmin, vmax=cs.vmax)
    #     # which does not work any more
    #     sm = plt.cm.ScalarMappable(norm=norm, cmap=self.color)
    #     sm.set_array([])
    #     self.fig.colorbar(sm, ax=self.ax, **self.colorbar_args)
    #     plt.subplots_adjust(wspace=0.05)
        # center the plot in the figure

    def is_3d(self) -> bool:
        return False

    def get_plot(self):
        return [self.ax.contourf(self.X, self.Y, self.Z, **self.plot_args)]

    def get_plot_args(self) -> dict[str, Any]:
        plt.style.use('_mpl-gallery-nogrid')
        self.num_levels = 20
        self.levels = np.linspace(self.Z_MIN, self.Z_MAX, self.num_levels)
        return {'cmap': self.color, 'levels': self.levels}

    def get_plot_args_extra(self):
        pass


class SurfaceContour(Plot):
    def __init__(self):
        super().__init__()
        self.offset_rate = 2
        self.plot_setup()
        # self.plot.pop()
        # self.num_plots -= 1
        self.ax.set_zlim(self.offset_rate*self.Z_MIN, self.Z_MAX)
        # self.custom_colorbar()
        self.show_plot()

    # def custom_colorbar(self):
    #     # remove the previous colorbar
    #     self.fig.delaxes(self.fig.axes[1])

    #     # create a new colorbar with self.levels
    #     norm = Normalize(vmin=self.Z_MIN, vmax=self.Z_MAX)
    #     # a previous version of this used
    #     # norm= matplotlib.colors.Normalize(vmin=cs.vmin, vmax=cs.vmax)
    #     # which does not work any more
    #     sm = plt.cm.ScalarMappable(norm=norm, cmap=self.color)
    #     sm.set_array([])
    #     self.fig.colorbar(sm, ax=self.ax, **self.colorbar_args)

    def is_3d(self) -> bool:
        return True

    def get_plot(self):
        # assume self.ax is Axes3D
        self.ax = cast(Axes3D, self.ax)
        return [self.ax.contour(self.X, self.Y, self.Z, **self.plot_args_extra), self.ax.plot_surface(self.X, self.Y, self.Z, **self.plot_args)]

    def get_plot_args(self) -> dict[str, Any]:
        # return {'rstride': 1, 'cstride': 1, 'facecolor': self.color_extra, 'linewidth': 0.01,
        #         'antialiased': True, 'color': 'w', 'shade': True, 'vmin': self.Z_MIN, 'vmax': self.Z_MAX, 'alpha': 0.9}
        return {'rstride': 1, 'cstride': 1, 'cmap': self.color, 'linewidth': 0.01,
                'antialiased': True, 'color': 'w', 'shade': True, 'vmin': self.Z_MIN, 'vmax': self.Z_MAX}

    def get_plot_args_extra(self):
        plt.style.use('_mpl-gallery-nogrid')
        self.num_levels = 20
        self.levels = np.linspace(self.Z_MIN, self.Z_MAX, self.num_levels)
        return {'zdir': 'z', 'offset': self.offset_rate * self.Z_MIN, 'cmap': self.color, 'levels': self.levels}


# count time
UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
# nu1 = float(sys.argv[1])
# nu2 = float(sys.argv[2])
# arg = int(sys.argv[3])
# Surface()
# elif arg == 2:
#     Contour(nu1, nu2)
# else:
SurfaceContour()
