"""
    Authors: Artur Przybysz (s202384), Elena Mongelli(s181214), Ivan Knezevic (s202386)
"""

# Exercise c

import math
import numpy as np
import matplotlib.pyplot as plt


# U(x) = e^x; 
# U_hat(x) = u_hat1 * N1 + u_hat2 * N2 + u_hat3 * N3

def u_x(x):
    return math.exp(x)


def N1(x, h1):
    if 0 <= x <= h1:
        return (h1 - x) / h1
    else:
        return 0


def N2(x, h1, h2):
    if 0 <= x <= h1:
        return x / h1
    if h1 <= x <= h1 + h2:
        return (h1 + h2 - x) / h2
    else:
        return 0


def N3(x, h1, h2):
    if h1 <= x <= h1 + h2:
        return (x - h1) / h2
    else:
        return 0


def u_hat(x, h1, h2, u_hat1, u_hat2, u_hat3):
    return u_hat1 * N1(x, h1) + u_hat2 * N2(x, h1, h2) + u_hat3 * N3(x, h1, h2)


def u_i(x, h1, h2):
    return u_x(0) * N1(x, h1) + u_x(h1) * N2(x, h1, h2) + u_x(h1 + h2) * N3(x, h1, h2)


if __name__ == '__main__':
    u_hat1 = 1
    u_hat2 = 2.6216
    u_hat3 = math.exp(2)
    h1 = 1
    h2 = 1

    u_hat_values = []
    u_x_values = []
    u_i_values = []
    x = np.arange(0, 2.1, 0.1)
    for i in x:
        u_hat_values.append(u_hat(i, h1, h2, u_hat1, u_hat2, u_hat3))
        u_x_values.append(u_x(i))
        u_i_values.append(u_i(i, h1, h2))

    plt.plot(x, u_hat_values, label="u_hat")
    plt.plot(x, u_x_values, label="u_x")
    plt.plot(x, u_i_values, label="u_i")
    plt.legend()
    plt.show()

    # Exercise d: redefing h1, h2 and u_hat2 since all the other values don't change

    h1 = 4 / 3
    h2 = 2 / 3
    u_hat2 = 3.6996

    u_hat_values = []
    u_x_values = []
    u_i_values = []
    x = np.arange(0, 2.1, 0.1)
    for i in x:
        u_hat_values.append(u_hat(i, h1, h2, u_hat1, u_hat2, u_hat3))
        u_x_values.append(u_x(i))
        u_i_values.append(u_i(i, h1, h2))

    plt.plot(x, u_hat_values, label="u_hat")
    plt.plot(x, u_x_values, label="u_x")
    plt.plot(x, u_i_values, label="u_i")
    plt.legend()
    plt.show()
