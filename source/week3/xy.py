import numpy as np


def xy(x0, y0, L1, L2, noelms1, noelms2, as_list=False):
    all_nodes = (1 + noelms2) * (1 + noelms1)
    X = np.zeros(all_nodes)

    for i, x in enumerate(np.linspace(start=x0, stop=x0 + L1, num=noelms2 + 1)):
        start_idx = i * (noelms1 + 1)
        stop_idx = start_idx + (noelms1 + 1)
        X[start_idx:stop_idx] = x

    first_row = np.linspace(start=y0, stop=y0 + L2, num=noelms1 + 1)[::-1]
    Y = np.tile(first_row, noelms2 + 1)
    if as_list:
        return X.tolist(), Y.tolist()
    return X, Y
