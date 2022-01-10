"""
    Authors: Artur Przybysz (s202384), Elena Mongelli(s181214), Ivan Knezevic (s202386)
"""
# 1.5d)

import numpy as np
import matplotlib.pyplot as plt

psi = 1
epsilon = [1, 0.01, 0.001]
x = np.arange(0, 1, 0.001)


if __name__ == '__main__':
    for ep in epsilon:
        # u=1/psi *((1+(np.exp(psi/ep)-1)*x-np.exp(x*psi/ep))/(np.exp(psi/ep)-1))

        a = float(psi / ep)
        # alternative form
        # u=1/psi*(1/(-1+np.exp(a)-1)-np.exp(x*a)/(-1+np.exp(a))+x)

        u = 1 / psi * ((np.exp(-a) + (x - x * np.exp(-a)) - np.exp(psi * (x - 1) / ep)) / (1 - np.exp(-a)))
        plt.plot(x, u)
        plt.show()
