import numpy as np

from source.week2.ex2_1 import xy, construct_element_table
from source.week2.ex2_2a import basfun

from scipy.sparse import spdiags, diags


def compute_k_rsn(abc, delta, r, s, lam1, lam2):
    br = abc[r, 1]
    bs = abc[s, 1]
    cr = abc[r, 2]
    cs = abc[s, 2]
    k_rsn = 1 / (delta * 4) * (lam1 * br * bs + lam2 * cr * cs)
    return k_rsn


def assembly(VX, VY, EToV, lam1, lam2, qt, M):
    """
        TODO
    """
    N = len(EToV)
    A = np.zeros((M, M))
    b = np.zeros(M)
    for n in range(1, N + 1):
        i, j, k = EToV[n]
        i -= 1
        j -= 1
        k -= 1

        q = qt[n]
        abc, delta = basfun(n, VX, VY, EToV)

        for r in range(3):
            b[i] += q[r]
            for s in range(3):
                k_rsn = compute_k_rsn(abc, delta, r, s, lam1, lam2)
                A[i, j] += k_rsn
    return A, b  # TODO: make sure A is sparse


def construct_qt(etov_dict, VX, VY):
    """
    Computes points as means of each elements' vertices.
    """
    qt_dict = dict()
    for e, (v1, v2, v3) in etov_dict.items():
        x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
        y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

        mean_1 = np.mean([x1, y1])
        mean_2 = np.mean([x2, y2])
        mean_3 = np.mean([x3, y3])

        qt_dict[e] = (mean_1, mean_2, mean_3)
    return qt_dict


def main():
    X, Y = xy(x0=0, y0=0, L1=1, L2=1, noelms1=4, noelms2=4)
    etov_dict, M = construct_element_table(4, 4)
    # TODO: get lambda value
    lam1 = 1
    lam2 = 1
    qt = construct_qt(etov_dict, X, Y)
    A, b = assembly(X, Y, etov_dict, lam1=lam1, lam2=lam2, qt=qt, M=M)
    # B, d = spdiags(A)
    print()


if __name__ == '__main__':
    main()
