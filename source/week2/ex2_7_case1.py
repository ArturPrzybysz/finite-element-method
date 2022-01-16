from collections import defaultdict

import numpy as np

from source.week2.ex2_1 import construct_element_table, xy
from source.week2.ex2_2b import outernormal
from source.week2.ex2_3 import assembly
from source.week2.ex2_6 import construct_boundary_edges, neumann_boundary_conditions
from source.week2.ex2_7 import dirichlet_boundary_conditions, find_boundary_indices, find_neighbors


def ut(x, y):
    return 3 * x + 5 * y - 7


def q(x, y, n1, n2):
    return -(lam1 * x * n1 + lam2 * y * n2)


def neumann_boundary_conditions_case1(VX, VY, EToV, gamma1, b1):
    b = np.array(b1)
    for e, r in gamma1:
        s = (r + 1) % 3
        vertice_tuple = EToV[e]
        global_r = vertice_tuple[r]
        global_s = vertice_tuple[s]
        x_r, y_r = VX[global_r - 1], VY[global_r - 1]
        x_s, y_s = VX[global_s - 1], VY[global_s - 1]

        xc = 0.5 * (x_r + x_s)
        yc = 0.5 * (y_r + y_s)
        n1, n2 = outernormal(e, r + 1, VX, VY, EToV)
        q12 = (q(xc, yc, n1, n2) / 2) * np.sqrt((x_r - x_s) ** 2 + (y_r - y_s) ** 2)
        b[global_r - 1] -= q12
        b[global_s - 1] -= q12
    return b


def dirichlet_boundary_conditions_case1(A_old, b_old, gamma2, EToV, VX, VY):
    A = np.array(A_old)
    b = np.array(b_old)
    S_gamma = find_boundary_indices(gamma2, EToV)

    for i in S_gamma:
        x_i, y_i = VX[i - 1], VY[i - 1]
        f_i = 3 * x_i + 5 * y_i - 7

        A[i - 1, i - 1] = 1
        b[i - 1] = f_i
        neighbors = find_neighbors(i, EToV)

        for j in neighbors:
            if j < i:
                b[j - 1] -= A[j - 1, i - 1] * f_i
                A[j - 1, i - 1] = 0
            else:
                b[j - 1] -= A[i - 1, j - 1] * f_i
                A[i - 1, j - 1] = 0
    return A, b


if __name__ == '__main__':
    noelms1 = 4
    noelms2 = 3
    L1 = 7.6
    L2 = 5.9
    x0 = -2.5
    y0 = -4.8
    lam1 = 1
    lam2 = 1
    tol = 5e-2
    # Solving the partial derivative of u_t
    qx = 3
    qy = 5
    qt = defaultdict(lambda: 0)

    # Building VX, VY and EToV based on the input information above. These % codes are explained in previous exercises. [VX, VY] = xy(x0, y0, L1, L2, noelms1, noelms2);
    VX, VY = xy(x0, y0, L1, L2, noelms1, noelms2)
    EToV, M = construct_element_table(noelms1, noelms2)
    A, b = assembly(VX, VY, EToV, lam1, lam2, qt, M)

    gamma2, gamma1 = construct_boundary_edges(VX, VY, EToV, tol=0.0005, x0=x0, y0=y0, L1=L1, L2=L2)
    b_nc = neumann_boundary_conditions_case1(VX, VY, EToV, gamma1, b)

    # Computing solution for Dirichlet boundary condition
    f = ut(VX, VY)
    # A_old, b_old, lower_edge, etov, VX, VY, exercise, test_case
    A_b, b_b = dirichlet_boundary_conditions_case1(A, b_nc, gamma2, EToV, VX, VY)
    u = np.linalg.solve(A_b, b_b)
    compare = np.stack([u, f])
    correct = np.isclose(u, f)

    print()
    # % Calculation the final solution and reshape it to 2-D array format
    # u = A_b \ b_b;
    # U = reshape(u,noelms2+1,noelms1+1);
    # disp("Solution in 2-D"); U
