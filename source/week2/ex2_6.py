"""
a) Write a Matlab routine that constructs an array for holding information about boundary
edges (beds) of the triangular mesh. (HINT: the use of a distance function fd can make the
job of finding nodes that belong to a boundary implicitly defined by the distance function
very easy. Below it is included in the argument list as an optional argument. top is a
tolerance to be used with the signed distance function to tune the selection of nodes.)
"""
import numpy as np
import matplotlib.pyplot as plt

from source.week2.ex2_1 import construct_element_table, xy
from source.week2.ex2_2a import basfun
from source.week2.ex2_2b import outernormal
from source.week2.ex2_3 import assembly

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


def q(VX, VY, etov, lam1, lam2, n, k, test_case=None):
    if test_case == 1:
        u_x = 3
        u_y = 5
        n = outernormal(n, k, VX, VY, etov)
        n1, n2 = n[0], n[1]
        return - lam1 * u_x * n1 - lam2 * u_y * n2

    if test_case == 2:
        u_xx = -np.sin(?) * np.sin(?)
        u_yy = -np.sin(?) * np.sin(?)
        return -u_xx - u_yy


def neumann_boundary_conditions(VX, VY, lam1, lam2, EToV, boundary_edges, qt, b1, test_case):
    b = np.array(b1)
    for e, r in boundary_edges:
        s = (r + 1) % 3
        vertice_tuple = EToV[e]
        # _, delta = basfun(e, VX, VY, EToV)
        global_r = vertice_tuple[r]
        global_s = vertice_tuple[s]
        x_r, y_r = VX[global_r - 1], VY[global_r - 1]
        x_s, y_s = VX[global_s - 1], VY[global_s - 1]

        q_1TODO = q(VX, VY, EToV, lam1, lam2, e, k=r, test_case=test_case)
        q_2TODO = q(VX, VY, EToV, lam1, lam2, e, k=s, test_case=test_case)
        q_1 = q_ / 2 * np.sqrt((x_r - x_s) ** 2 + (y_r - y_s) ** 2)  # q_2 == q_1
        b[global_r - 1] -= q_1
        b[global_s - 1] -= q_1
    return b


def construct_qt_2_7(etov, VX, VY, test_case, lam1=1, lam2=1):
    qt_dict = dict()
    for e, (v1, v2, v3) in etov.items():
        if test_case == 1:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            # def f(x, y, lam1, lam2, n, k):
            #     u_x = 3
            #     u_y = 5
            #     n = outernormal(n, k, VX, VY, etov)
            #     n1, n2 = n[0], n[1]
            #     return - lam1 * u_x * n1 - lam2 * u_y * n2

            def f(x, y, lam1, lam2, n, k):
                u_xx = 0
                u_yy = 0
                return - u_xx - u_yy

            q1 = f(x1, y1, lam1=lam1, lam2=lam2, n=e, k=1)
            q2 = f(x2, y2, lam1=lam1, lam2=lam2, n=e, k=2)
            q3 = f(x3, y3, lam1=lam1, lam2=lam2, n=e, k=3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
        elif test_case == 2:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            def f(x, y, lam1, lam2, n, k):
                u_xx = -np.sin(x) * np.sin(y)
                u_yy = -np.sin(x) * np.sin(y)
                return -u_xx - u_yy

            # def f(x, y, lam1, lam2, n, k):
            #     u_x = np.cos(x) * np.sin(y)
            #     u_y = np.sin(x) * np.cos(y)
            #     n = outernormal(n, k, VX, VY, etov)
            #     n1, n2 = n[0], n[1]
            #     return - lam1 * u_x * n1 - lam2 * u_y * n2

            q1 = f(x1, y1, lam1=lam1, lam2=lam2, n=e, k=1)
            q2 = f(x2, y2, lam1=lam1, lam2=lam2, n=e, k=2)
            q3 = f(x3, y3, lam1=lam1, lam2=lam2, n=e, k=3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
    return qt_dict


def test_case_data_ex_2_7(test_case):
    if test_case == 1:
        x0 = -2.5
        y0 = -4.8
        L1 = 7.6
        L2 = 5.9
        X, Y = xy(x0=x0, y0=y0, L1=L1, L2=L2, noelms1=4, noelms2=3)
        etov, M = construct_element_table(4, 3)
        lam1 = 1
        lam2 = 1
        qt = construct_qt_2_7(etov, X, Y, test_case=test_case)
        return test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2
    elif test_case == 2:
        x0 = -2.5
        y0 = -4.8
        L1 = 7.6
        L2 = 5.9
        noelms1 = 32
        noelms2 = 32
        X, Y = xy(x0=x0, y0=y0, L1=L1, L2=L2, noelms1=noelms1, noelms2=noelms2)
        etov, M = construct_element_table(noelms1, noelms2)
        lam1 = 1
        lam2 = 1
        qt = construct_qt_2_7(etov, X, Y, test_case=test_case)
        return test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2
    else:
        raise Exception(f"Unknown test case: {test_case}")


def main():
    nr_of_test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2 = test_case_data_ex_2_7(1)
    lower_points, upper_points = construct_boundary_edges(X, Y, etov, tol=0.0005, x0=x0, y0=y0, L1=L1, L2=L2)
    visualise(X, Y, etov, lower_points)
    visualise(X, Y, etov, upper_points)
    A, b = assembly(X, Y, etov, lam1=lam1, lam2=lam2, qt=qt, M=M)
    b_neumann1 = neumann_boundary_conditions(X, Y, lam1, lam2, etov, upper_points, qt, b, test_case=nr_of_test_case)
    b_neumann2 = neumann_boundary_conditions(X, Y, lam1, lam2, etov, lower_points, qt, b_neumann1,
                                             test_case=nr_of_test_case)
    u = np.linalg.solve(A, b)
    print()



if __name__ == '__main__':
    main()
