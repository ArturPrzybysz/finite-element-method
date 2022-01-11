import numpy as np

from source.week2.ex2_1 import xy, construct_element_table
from source.week2.ex2_2a import basfun

from scipy.sparse import spdiags, diags


def compute_k_rsn(abc, delta, r, s, lam1, lam2):
    br = abc[r, 1]
    bs = abc[s, 1]
    cr = abc[r, 2]
    cs = abc[s, 2]
    k_rsn = (lam1 * br * bs + lam2 * cr * cs) / np.abs(delta * 4)
    return k_rsn


def assembly(VX, VY, EToV, lam1, lam2, qt, M):
    """
        TODO: explanation
    """
    N = len(EToV)
    iL = np.zeros(3 * 3 * N)
    jL = np.zeros(3 * 3 * N)
    sL = np.zeros(3 * 3 * N)
    count = 0
    A = np.zeros((M, M))
    b = np.zeros(M)
    for n in range(1, N + 1):
        # q = qt[n]
        abc, delta = basfun(n, VX, VY, EToV)
        for r in range(3):
            i = EToV[n][r] - 1
            qr = qt[n]
            b[i] += qr
            for s in range(3):
                j = EToV[n][s] - 1
                k_rsn = compute_k_rsn(abc, delta, r, s, lam1, lam2)
                iL[count] = i
                jL[count] = j
                sL[count] = k_rsn
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


def main():
    X, Y = xy(x0=0, y0=0, L1=1, L2=1, noelms1=4, noelms2=3)
    etov_dict, M = construct_element_table(4, 3)
    # TODO: get lambda value
    lam1 = 1
    lam2 = 1
    qt = construct_qt(etov_dict, X, Y, test_case=1)
    A, b = assembly(X, Y, etov_dict, lam1=lam1, lam2=lam2, qt=qt, M=M)
    # B, d = spdiags(A)
    print()


if __name__ == '__main__':
    main()
