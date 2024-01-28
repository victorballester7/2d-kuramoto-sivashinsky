import numpy as np
from read_data import read_data_energy_return
import time
import sys

UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
filename_E = "data/energy_return.txt"
data_E = read_data_energy_return(filename_E)

try:
    t_min = float(sys.argv[1])
except IndexError:
    t_min = 200.0

try:
    idx_1 = np.where(data_E[:, 0] > t_min)[0][0]
except IndexError:
    idx_1 = 0
times = data_E[idx_1:, 0]
data_E = data_E[idx_1:, 1]

E_last = data_E[-1]
data_E = np.abs(data_E - E_last)

EPS = 0.5
data_E = (data_E < EPS)

# get the indices where data_E is True
idx = np.where(data_E)[0]
times = times[idx]

for i in range(len(times) - 1):
    print(f"Estimated period {i + 1}: {times[i + 1] - times[i]}")

end_time = time.time()
print("Total time for estimating the period: ", int(
    UNIT_TIME * (end_time - start_time)), LABEL_TIME)
