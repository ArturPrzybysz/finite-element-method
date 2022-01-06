# 1.5d)

import numpy as np
import matplotlib.pyplot as plt

psi = 1
epsilon = [1, 0.01, 0.001]
x = np.arange(0, 1, 0.001)
if __name__ == '__main__':

    for ep in epsilon:
        # u=1/psi *((1+(np.exp(psi/ep)-1)*x-np.exp(x*psi/ep))/(np.exp(psi/ep)-1))
        # alternative form
        u = 1 / psi * (1 / (-1 + np.exp(psi / ep) - 1) - np.exp(x * psi / ep) / (-1 + np.exp(psi / ep)) + x)
        # print(u)
        plt.plot(x, u)
        # TODO Title etc.
        plt.show()
