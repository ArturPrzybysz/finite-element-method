import numpy as np


def sparse_diags(A):
    # x = dia_matrix(A)
    max_diag_size = A.diagonal(0).shape[0]
    d = []
    data = []
    for diag in range(-A.shape[0], A.shape[1]):
        diag_value = A.diagonal(diag)
        if np.any(diag_value):
            d.append(diag)
            if diag < 0:
                x = np.zeros(max_diag_size)
                x[np.abs(diag):] = diag_value
                diag_value = x

            if diag > 0:
                x = np.zeros(max_diag_size)
                x[:-diag] = diag_value
                diag_value = x

            data.append(diag_value)

    d = np.array(d)
    B = np.array(data)
    return B, d


if __name__ == '__main__':
    from scipy.sparse import dia_matrix

    n = 10
    ex = np.ones(n)
    data = np.array([ex, 2 * ex, ex])
    offsets = np.array([-1, 0, 1])
    A = dia_matrix((data, offsets), shape=(n, n)).toarray()
    data2, offsets2 = sparse_diags(A)
    A2 = dia_matrix((data2, offsets2), shape=(n, n)).toarray()
    assert np.allclose(A, A2)
