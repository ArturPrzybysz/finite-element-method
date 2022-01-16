import numpy as np

from source.week2a.basfun import basfun


def compute_k_rsn(abc, delta, r, s, lam1, lam2):
    br = abc[r, 1]
    bs = abc[s, 1]
    cr = abc[r, 2]
    cs = abc[s, 2]
    k_rsn = (lam1 * br * bs + lam2 * cr * cs) / np.abs(delta * 4)
    return k_rsn


def assembly(VX, VY, EToV, lam1, lam2, qt, M):
    N = len(EToV)
    A = np.zeros((M, M))
    b = np.zeros(M)
    for n in range(1, N + 1):
        abc, delta = basfun(n, VX, VY, EToV)
        for r in range(3):
            i = EToV[n][r] - 1
            qr = qt(n) * np.abs(delta) / 3  # (2.28)
            b[i] += qr
            for s in range(3):
                j = EToV[n][s] - 1
                k_rsn = compute_k_rsn(abc, delta, r, s, lam1, lam2)  # (2.27)
                A[i, j] += k_rsn
    return A, b
