import numpy as np

from source.week2.ex2_1 import xy, construct_element_table
from source.week2.ex2_2a import basfun

from scipy.sparse import dia_matrix


def compute_k_rsn(abc, delta, r, s, lam1, lam2):
    br = abc[r, 1]
    bs = abc[s, 1]
    cr = abc[r, 2]
    cs = abc[s, 2]
    k_rsn = (lam1 * br * bs + lam2 * cr * cs) / np.abs(delta * 4)
    return k_rsn


def assembly(VX, VY, EToV, lam1, lam2, qt, M, ):
    """
        TODO: explanation
    """
    N = len(EToV)
    A = np.zeros((M, M))
    b = np.zeros(M)
    for n in range(1, N + 1):
        abc, delta = basfun(n, VX, VY, EToV)
        for r in range(3):
            i = EToV[n][r] - 1

            qr = qt[n] * np.abs(delta) / 3
            b[i] += qr
            for s in range(3):
                j = EToV[n][s] - 1
                k_rsn = compute_k_rsn(abc, delta, r, s, lam1, lam2)
                A[i, j] += k_rsn
    return A, b  # TODO: make sure A is sparse


def construct_qt(etov_dict, VX, VY, test_case=None):
    """
    Computes points as means of each elements' vertices.
    """
    qt_dict = dict()
    for e, (v1, v2, v3) in etov_dict.items():
        if test_case == 1:
            qt_dict[e] = 0
        elif test_case == 2:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            def f(x, y):
                return -6 * x + 2 * y - 2

            q1 = f(x1, y1)
            q2 = f(x2, y2)
            q3 = f(x3, y3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
        else:
            raise Exception("Unknown test case")
    return qt_dict


def test_case_data_ex_2_3(nr_of_test_case):
    if nr_of_test_case == 1:
        x0 = 0
        y0 = 0
        L1 = 1
        L2 = 1
        X, Y = xy(x0=x0, y0=y0, L1=L1, L2=L2, noelms1=4, noelms2=3)
        etov, M = construct_element_table(4, 3)
        # TODO: get lambda value
        lam1 = 1
        lam2 = 1
        qt = construct_qt(etov, X, Y, test_case=1)
        return nr_of_test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2
    elif nr_of_test_case == 2:
        x0 = -2.5
        y0 = -4.8
        L1 = 7.6
        L2 = 5.9
        X, Y = xy(x0=x0, y0=y0, L1=L1, L2=L2, noelms1=4, noelms2=3)
        etov, M = construct_element_table(4, 3)
        lam1 = 1
        lam2 = 1
        qt = construct_qt(etov, X, Y, test_case=2)
        return nr_of_test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2
    else:
        raise Exception("Unknown test case")


def sparse_diags(A):
    # x = dia_matrix(A)
    max_diag_size = A.diagonal(0).shape[0]
    d = []
    data = []
    for diag in range(-A.shape[0], A.shape[1]):
        diag_value = A.diagonal(diag)
        if np.any(diag_value):
            d.append(diag)

            if diag > 0:
                x = np.zeros(max_diag_size)
                x[diag:] = diag_value
                diag_value = x
            data.append(diag_value)

    d = np.array(d)
    B = np.zeros((len(d), max_diag_size))
    return B, d


def main():
    nr_of_test_case, X, Y, etov_dict, M, lam1, lam2, qt, x0, y0, L1, L2 = test_case_data_ex_2_3(2)
    A, b = assembly(X, Y, etov_dict, lam1=lam1, lam2=lam2, qt=qt, M=M)
    B, d = sparse_diags(A)

    print()


if __name__ == '__main__':
    main()
