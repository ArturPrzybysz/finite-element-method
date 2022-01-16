import numpy as np


def basfun(element_idx, X, Y, etov_dict):
    """
    :param element_idx:
    :param X:
    :param Y:
    :param etov_dict:
    :return:
    """

    # Formula 2.3

    abc = np.zeros((3, 3))
    v1, v2, v3 = etov_dict[element_idx]
    v1 -= 1
    v2 -= 1
    v3 -= 1

    abc[0][0] = X[v2] * Y[v3] - X[v3] * Y[v2]
    #   1  b
    abc[0][1] = Y[v2] - Y[v3]
    #   1  c
    abc[0][2] = X[v3] - X[v2]

    #   2  a
    abc[1][0] = X[v1] * Y[v3] - X[v3] * Y[v1]
    #   2  b
    abc[1][1] = -(Y[v1] - Y[v3])
    #   2  c
    abc[1][2] = X[v3] - X[v1]

    #   3  a
    abc[2][0] = X[v1] * Y[v2] - X[v2] * Y[v1]
    #   3  b
    abc[2][1] = Y[v1] - Y[v2]
    #   3  c
    abc[2][2] = X[v2] - X[v1]

    # Formula 2.1 TODO
    # Delta : positive due to clockwise ordering
    _delta = 0.5 * (X[v2] * Y[v3] - Y[v2] * X[v3] - (X[v1] * Y[v3] - Y[v1] * X[v3]) + X[v1] * Y[v2] - Y[v1] * X[v2])

    return abc, _delta
