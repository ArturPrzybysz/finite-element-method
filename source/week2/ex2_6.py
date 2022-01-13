"""
a) Write a Matlab routine that constructs an array for holding information about boundary
edges (beds) of the triangular mesh. (HINT: the use of a distance function fd can make the
job of finding nodes that belong to a boundary implicitly defined by the distance function
very easy. Below it is included in the argument list as an optional argument. top is a
tolerance to be used with the signed distance function to tune the selection of nodes.)
"""
import numpy as np
import matplotlib.pyplot as plt

from source.week2.ex2_2a import basfun
from source.week2.ex2_3 import test_case_data_ex_2_3, construct_qt, assembly

"""
Hint: Check slide 30 in lecture 4 on Mesh Generation. Signed distance functions are available in the DistMesh package
"""


def visualise(VX, VY, EToV, boundary_edges):
    for e, indices in EToV.items():
        for r, v in enumerate(indices):
            x, y = VX[v - 1], VY[v - 1]
            plt.scatter(x, y, c='grey')
    for idx in range(len(boundary_edges)):
        e, r = boundary_edges[idx]
        v = EToV[e][r]
        x, y = VX[v - 1], VY[v - 1]
        plt.scatter(x, y, c='red')
    plt.show()


def construct_boundary_edges(VX, VY, EToV, tol, x0, y0, L1, L2):
    upper_points = set()
    lower_points = set()
    for e, indices in EToV.items():
        for r, v in enumerate(indices):
            x, y = VX[v - 1], VY[v - 1]
            if np.allclose(x, x0, atol=tol) or np.allclose(y, y0 + L2, atol=tol):
                lower_points.add((e, r))
            if np.allclose(y, y0, atol=tol) or np.allclose(x, x0 + L1, atol=tol):
                upper_points.add((e, r))

    # for idx, overlap in enumerate(upper_points.intersection(lower_points)):
    #     if idx % 2 == 0:
    #         upper_points.remove(overlap)
    #     else:
    #         lower_points.remove(overlap)
    return list(lower_points), list(upper_points)


def neumann_boundary_conditions(VX, VY, EToV, boundary_edges, qt, b):
    for e, r in boundary_edges:
        s = (r + 1) % 3
        vertice_tuple = EToV[e]
        _, delta = basfun(e, VX, VY, EToV)
        global_r = vertice_tuple[r]
        global_s = vertice_tuple[s]
        x_r, y_r = VX[global_r - 1], VY[global_r - 1]
        x_s, y_s = VX[global_s - 1], VY[global_s - 1]
        q = qt[e] * np.abs(delta) / 3
        q_1 = q / 2 * np.sqrt((x_r - x_s) ** 2 + (y_r - y_s) ** 2)  # q_2 == q_1
        b[global_r - 1] -= q_1
        b[global_s - 1] -= q_1
    return b


def main():
    nr_of_test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2 = test_case_data_ex_2_3(2)
    lower_points, upper_points = construct_boundary_edges(X, Y, etov, tol=0.0005, x0=x0, y0=y0, L1=L1, L2=L2)
    visualise(X, Y, etov, lower_points)
    visualise(X, Y, etov, upper_points)
    A, b = assembly(X, Y, etov, lam1=lam1, lam2=lam2, qt=qt, M=M)
    b_neumann1 = neumann_boundary_conditions(X, Y, etov, upper_points, qt, b)
    # b_neumann2 = neumann_boundary_conditions(X, Y, etov, lower_points, qt, b_neumann1)
    u = np.linalg.solve(A, b)
    print()


if __name__ == '__main__':
    main()
