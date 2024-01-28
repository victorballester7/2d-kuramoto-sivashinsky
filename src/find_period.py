import numpy as np
from read_data import read_data_energy
import time
import sys

UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
filename_E = "data/energy.txt"
data_u = read_data_energy(filename_E)

try:
    t_min = float(sys.argv[1])
except IndexError:
    t_min = 200.0

try:
    idx_1 = np.where(data_u[:, 0] > t_min)[0][0]
except IndexError:
    idx_1 = 0
times = data_u[idx_1:, 0]
data_u = data_u[idx_1:, 3]

U_MAX = np.max(data_u)
U_MIN = np.min(data_u)

L = (U_MAX + U_MIN) / 2
# (we substract an average value to the data_u to ensure there will be 0 in the image of data_u)
data_u = data_u - L

# data_u will contain points of a sinusoidal function. We want to find the
# period of this function.
data_u_shifted = np.roll(data_u, 1)
data_u_shifted[0] = data_u_shifted[1]

# We want to find the points where the sign of data_u changes from + to -.
change_sign = (np.sign(data_u) != np.sign(
    data_u_shifted)) * (np.sign(data_u_shifted) == 1)

# we want to store the values of times where u is = 0. This is done using
# 1st order linear interpolation. And we get the time t = [(y1 * t0 - y0 *
# t1)] / (y1 - y0)
idx = np.where(change_sign)[0]
times = ((data_u[idx] * times[idx - 1] - data_u[idx - 1]
         * times[idx])) / (data_u[idx] - data_u[idx - 1])

periods = np.diff(times)

length = len(periods)
periods_shown = 10
if length < periods_shown:
    periods_shown = length
for i in range(periods_shown):
    aux = - periods_shown + i
    print(f"Estimated period {length +aux + 1}: {periods[aux]}")

end_time = time.time()
print("Total time for estimating the period: ", int(
    UNIT_TIME * (end_time - start_time)), LABEL_TIME)
