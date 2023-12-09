import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from matplotlib import colormaps as cm
from matplotlib.ticker import LinearLocator
import os
from mpl_toolkits.mplot3d import Axes3D
from typing import cast


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
  ax = cast(Axes3D, fig.add_subplot(projection='3d'))

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
  # Z = np.sin(2*X) + np.cos(2*Y+headers[0])

  # DEFAULTS
  FONTSIZE_TITLE = 15
  FONTSIZE_TIME = 10

  # PLOT TYPE
  # Use the appropriate number for each plot type
  # 1 - surface
  # 2 - contour
  PLOT_TYPE = 1 
  color = cm['viridis']
  FPS = 50
  Z_MAX = np.max(np.abs(DATA_PLOT))
  Z_MIN = -Z_MAX
  plot_args = {}
  time_text = ax.text2D(0.5, 0.95, '', transform=ax.transAxes,fontsize=FONTSIZE_TIME, horizontalalignment='center')
  
  # ax.zaxis.set_major_locator(LinearLocator())
  # A StrMethodFormatter is used automatically
  # ax.zaxis.set_major_formatter('{x:.02f}')

  # Customize axis
  ax.set_title('2D Kuramoto-Sivashinsky equation',fontsize=FONTSIZE_TITLE, fontweight='bold')
  # center text
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')
  ax.set_zlim(Z_MIN, Z_MAX)

  # Plot the surface.

  plot = []
  if PLOT_TYPE == 1:
    plot_args = {'rstride': 1, 'cstride': 1, 'cmap': color,'linewidth': 0.01, 'antialiased': True, 'color': 'w', 'shade': True, 'vmin': Z_MIN, 'vmax': Z_MAX}
    plot.append(ax.plot_surface(X, Y, Z, **plot_args))
  elif PLOT_TYPE == 2:
    plot_args = {'cmap': color,'vmin': Z_MIN, 'vmax': Z_MAX}
    plot.append(ax.imshow(Z, **plot_args))
  # surf = ax.plot_surface(X, Y, Z, **plot_args)
  # Add a color bar which maps values to colors.
  cb = fig.colorbar(plot[0], shrink=0.5, aspect=10,location='left')
  cax = cb.ax

  # Animation update function
  def init():
    return plot[0], time_text
  
  def update(frame):
    # Clear the previous frame
    # ax.clear()
    plot[0].remove()

    # Update the arrays
    Z = data_blocks[frame]
    # Z = np.sin(2*X) + np.cos(2*Y+10*headers[frame])

    # Update the time text
    time_text.set_text('t = %.2f' % headers[frame])

    # Update the scatter plot
    if PLOT_TYPE == 1:
      plot[0] = ax.plot_surface(X, Y, Z, **plot_args)
    elif PLOT_TYPE == 2:
      plot[0] = ax.imshow(Z, **plot_args)
    # update the colorbar
    # fig.colorbar(surf[0], cax=cax)
    # cb.update_normal(surf[0])

    return plot[0], time_text

  # # Create the animation
  ani = FuncAnimation(fig, update, frames=len(data_blocks), interval=1000/FPS, blit=False, init_func=init)
  # Show the animation
  plt.show()
  

plot()
