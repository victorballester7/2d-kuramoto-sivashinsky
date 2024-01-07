import matplotlib.pyplot as plt
import numpy as np


def plotFunction(argx, argy):
    if type(argx) == list or type(argx) == np.ndarray and type(argy) == list or type(argy) == np.ndarray:
        plt.plot(argx, argy)
        plt.grid()
        plt.show()
    else:
        print("Not all arguments passed are lists. Call failed.")
