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
from read_data import read_data_file, read_data_file_freq


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
    color = cm['inferno']
    color_extra = 'royalblue'
    colorbar_args = {'shrink': 0.5, 'aspect': 10, 'location': 'left'}
    FPS = 50
    interval = 1000 / FPS  # interval between frames in milliseconds
    duration = 40  # duration of the animation in seconds

    def __init__(self):
        pass

    def plot_setup(self, isfreq=False, t_min=0, save=False):
        # Create a figure
        self.fig = plt.figure()
        self.ax = self.get_ax()

        # Read data
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Replace with the actual file path
        file_path = script_dir + '/../data/solution.txt'
        file_path_freq = script_dir + '/../data/freq.txt'
        self.headers, self.data_blocks = read_data_file(file_path)
        self.headers_freq, self.data_blocks_freq = read_data_file_freq(
            file_path_freq)
        EPS = 0.000001
        try:
            idx_1 = np.where(self.headers > t_min - EPS)[0][0]
        except IndexError:
            idx_1 = 0
        try:
            idx_2 = np.where(self.headers_freq > t_min - EPS)[0][0]
        except IndexError:
            idx_2 = 0

        self.headers = self.headers[idx_1:]
        self.data_blocks = self.data_blocks[idx_1:]
        self.headers_freq = self.headers_freq[idx_2:]
        self.data_blocks_freq = self.data_blocks_freq[idx_2:]

        # in practice it takes more time to animate, so we divide the duration
        # by 2 to keep the seconds more or less the same
        if not save:
            self.duration /= 2

        # frames to be animated
        self.num_frames = int(self.duration * 1000.0 / self.interval)
        if isfreq:
            if self.num_frames > len(self.data_blocks_freq):
                self.num_frames = len(self.data_blocks_freq)
        else:
            if self.num_frames > len(self.data_blocks):
                self.num_frames = len(self.data_blocks)

        # Speed of the animation
        if isfreq:
            self.speed = len(self.data_blocks_freq) / self.num_frames
        else:
            self.speed = len(self.data_blocks) / self.num_frames

        print("Number of frames: ", self.num_frames)
        print("Speed: ", self.speed)
        print("Interval: ", self.interval)
        print("Expected duration: ", self.duration)
        # print("len(self.data_blocks): ", len(self.data_blocks))
        # print("len(self.headers): ", len(self.headers))
        # print("len(self.data_blocks_freq): ", len(self.data_blocks_freq))

        # Create data
        self.nx = self.data_blocks.shape[1]
        self.ny = self.data_blocks.shape[2]

        self.X = np.linspace(0, 2 * np.pi, self.nx)
        self.Y = np.linspace(0, 2 * np.pi, self.ny)
        self.X, self.Y = np.meshgrid(self.X, self.Y, indexing='ij')
        if isfreq:
            self.Z = self.data_blocks_freq[0]
            self.Z_MAX = np.max(self.data_blocks_freq)
            self.Z_MIN = 0
        else:
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
        # subtitle, fontsize=self.FONTSIZE_TITLE - 2, y=self.y_pos_title -
        # 0.05)

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
        if isfreq:
            self.fig.colorbar(self.plot[0], **self.colorbar_args)
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

            real_frame = int(frame * self.speed)

            # Update the arrays
            if isfreq:
                self.Z = self.data_blocks_freq[real_frame]
            else:
                self.Z = self.data_blocks[real_frame]
            # Z = np.sin(2*X) + np.cos(2*Y+10*headers[frame])

            # Update the time text
            self.time_text.set_text('t = %.2f' %
                                    self.headers[real_frame])

            # Update the plot
            # self.plot[0] = self.get_plot()
            tmp = self.get_plot()
            for i in range(self.num_plots):
                self.plot[i] = tmp[i]

            return self.plot[0], self.time_text

        # Create the animation
        self.ani = FuncAnimation(self.fig, update, frames=self.num_frames,
                                 interval=self.interval, blit=False, init_func=init)

        end_time = time.time()
        print("Total time for animating: ", int(
            UNIT_TIME * (end_time - start_time)), LABEL_TIME)

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
    def __init__(self, t_min, save):
        super().__init__()
        self.offset_rate = 2
        self.plot_setup(t_min=t_min, save=save)
        # self.plot.pop()
        # self.num_plots -= 1
        self.ax.set_zlim(self.offset_rate * self.Z_MIN, self.Z_MAX)
        # self.custom_colorbar()
        if save:
            # save all the frames of the animation
            filename = 'latex/presentation/videos/animation_' + \
                str(nu1) + '_' + str(nu2) + '.mp4'
            self.ani.save(filename, writer='ffmpeg', dpi=300)
        else:
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
        return [self.ax.contour(self.X, self.Y, self.Z, **self.plot_args_extra),
                self.ax.plot_surface(self.X, self.Y, self.Z, **self.plot_args)]

    def get_plot_args(self) -> dict[str, Any]:
        # return {'rstride': 1, 'cstride': 1, 'facecolor': self.color_extra, 'linewidth': 0.01,
        #         'antialiased': True, 'color': 'w', 'shade': True, 'vmin': self.Z_MIN, 'vmax': self.Z_MAX, 'alpha': 0.9}
        return {'rstride': 1, 'cstride': 1, 'cmap': self.color, 'linewidth': 0.01,
                'antialiased': True, 'color': 'w', 'shade': True, 'vmin': self.Z_MIN, 'vmax': self.Z_MAX}

    def get_plot_args_extra(self):
        plt.style.use('_mpl-gallery-nogrid')
        self.num_levels = 20
        self.levels = np.linspace(self.Z_MIN, self.Z_MAX, self.num_levels)
        return {'zdir': 'z', 'offset': self.offset_rate *
                self.Z_MIN, 'cmap': self.color, 'levels': self.levels}


class Heatmap(Plot):
    def __init__(self, t_min):
        super().__init__()
        self.plot_setup(isfreq=True, t_min=t_min)
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
        return [self.ax.imshow(self.Z, **self.plot_args)]

    def get_plot_args(self) -> dict[str, Any]:
        plt.style.use('_mpl-gallery-nogrid')
        return {'cmap': self.color}

    def get_plot_args_extra(self):
        pass


class FreqPlot(Plot):
    def __init__(self):
        super().__init__()
        self.plot_setup(True)
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
        return True

    def get_plot(self):
        ys = np.arange(len(self.Z[0]))
        return [self.ax.plot(ys, self.Z[i], zs=i, zdir='x',
                             **self.plot_args) for i in range(len(self.Z))]

    def get_plot_args(self) -> dict[str, Any]:
        # plt.style.use('_mpl-gallery-nogrid')
        return {'color': 'red'}

    def get_plot_args_extra(self):
        pass


# count time
UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
i = int(sys.argv[1])
try:
    t_min = float(sys.argv[2])
    nu1 = float(sys.argv[3])
    nu2 = float(sys.argv[4])
except IndexError:
    t_min = 0.0
    nu1 = np.nan
    nu2 = np.nan
if i == 1 or i == 3:
    SurfaceContour(t_min, True)  # save activated
# restart all plt config
plt.rcdefaults()
if i == 2 or i == 3:
    # FreqPlot()
    Heatmap(t_min)
