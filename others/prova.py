from __future__ import division
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import numpy as np

plot_args = {'rstride': 1, 'cstride': 1, 'cmap':
             cm.bwr, 'linewidth': 0.01, 'antialiased': True, 'color': 'w',
             'shade': True}

size = 100
soln = np.zeros((size, size))
midpoint = size // 2
soln[midpoint, midpoint] = 1

# first frame
X = range(size)
Y = range(size)
X, Y = np.meshgrid(X, Y)
plot = ax.plot_surface(X, Y, soln, **plot_args)
pam_ani = animation.FuncAnimation(fig, data_gen, fargs=(soln, plot),
                                  interval=30, blit=False)


def data_gen(framenumber, soln, plot):
    # change soln variable for the next frame
    ...
    ax.clear()
    plot = ax.plot_surface(X, Y, soln, **plot_args)
    return plot,


plt.show()
