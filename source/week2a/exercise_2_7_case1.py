from collections import defaultdict

import numpy as np

from source.week2a.assembly import assembly
from source.week2a.conelmtab import construct_element_table
from source.week2a.construct_beds import construct_boundary_edges
from source.week2a.dirbc import dirbc
from source.week2a.neubc import neumann_boundary_conditions
from source.week2a.xy import xy
import matplotlib.pyplot as plt

noelms1 = 4
noelms2 = 3
L1 = 7.6
L2 = 5.9
x0 = -2.5
y0 = -4.8
lam1 = 1
lam2 = 1
tol = 5e-2

ut = lambda x: 3 * x[0] + 5 * x[1] - 7

# Solving the partial derivative of u_t
qx = 3
qy = 5
qt = lambda _: 0


def visualise(VX, VY, EToV, boundary_edges, boundary_vertices):
    for e, indices in EToV.items():
        for r, v1 in enumerate(indices):
            x, y = VX[v1 - 1], VY[v1 - 1]
            plt.scatter(x, y, c='grey')
    for idx in range(len(boundary_edges)):
        e, r = boundary_edges[idx]
        s = (r + 1) % 3
        v1 = EToV[e][r]
        v2 = EToV[e][s]
        x_r, y_r = VX[v1 - 1], VY[v1 - 1]
        x_s, y_s = VX[v2 - 1], VY[v2 - 1]
        plt.plot([x_r, x_s], [y_r, y_s], c='green', alpha=0.3)
    for vertice in boundary_vertices:
        x_r, y_r = VX[vertice - 1], VY[vertice - 1]
        plt.scatter(x_r, y_r, c='black', alpha=1)

    plt.show()


if __name__ == '__main__':
    [VX, VY] = xy(x0, y0, L1, L2, noelms1, noelms2)
    [EToV, M] = construct_element_table(noelms1, noelms2)
    [A, b] = assembly(VX, VY, EToV, lam1, lam2, qt, M)
    # Computing the beds (the index of element and index of the face of boundary nodes)
    gamma1, gamma2, gnodes1, gnodes2 = construct_boundary_edges(VX, VY, EToV, tol, x0, y0, L1, L2)

    # % Creating the q function in the formula from exercise 2.6 for Neumann boundary condition
    # q = @(x,y,n1,n2) -(lam1*qx*n1 + lam2*qy*n2);

    def q(x, y, n1, n2):
        return -(lam1 * qx * n1 + lam2 * qy * n2)


    # % Computing solution for Neumann boundary condition
    # b_nc = neubc(VX,VY,EToV,Gamma1,q,b);
    b_nc_expected = [2.95, -5.90, -5.90, -7.70, 0, 0, 0,
                     -9.50, 0, 0, 0, -9.50, 0, 0, 0, -9.50, 0, 0, 0,
                     -4.75]
    b_nc = neumann_boundary_conditions(VX, VY, lam1, lam2, EToV, gamma1, q, b, b_nc_expected)

    b_nc_compare = np.stack([b_nc, b_nc_expected])
    b_nc_correct = b_nc == b_nc_expected
    # % Computing solution for Dirichlet boundary condition
    f = ut((VX, VY))
    Gamma1 = list(sorted([EToV[e][v] for e, v in gamma1]))
    Gamma2 = list(sorted([EToV[e][v] for e, v in gamma2]))
    gnodes1 = list(sorted(gnodes1))
    gnodes2 = list(sorted(gnodes2))
    [A_b, b_dir] = dirbc(gnodes1, f, A, b)

    b_dir_target = np.array(
        [-9.00, -10.2474576271186, -5.90, -7.70, -3.30, -3.18813559322034, 0, -9.50, 2.40, 2.31864406779661, 0, -9.50,
         8.10, 11.9312716820299, -6.07251461988304, -17.6254385964912, 13.8000000000000, 3.96666666666667,
         -5.86666666666667, -15.7000000000000])
    compare_b_b = np.stack([b_dir, b_dir_target])
    b_b_correct = np.isclose(b_dir, b_dir_target)
    print()
    u_hat = np.linalg.solve(A_b, b_dir)
    compare = np.stack([u_hat, f])
    correct = np.isclose(u_hat, f)
    print()
    #
    # % Calculation the final solution and reshape it to 2-D array format
    # u = A_b \ b_b;
    # U = reshape(u,noelms2+1,noelms1+1);
    # disp("Solution in 2-D"); U
