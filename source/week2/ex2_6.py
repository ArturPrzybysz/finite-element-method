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
from source.week2.ex2_4 import test_case_data

"""
Hint: Check slide 30 in lecture 4 on Mesh Generation. Signed distance functions are available in the DistMesh package
"""


def visualise(VX, VY, EToV, boundary_edges):
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
        plt.scatter(x_r, y_r, c='red')
        plt.plot([x_r, x_s], [y_r, y_s], c='green')

    plt.show()


def construct_boundary_edges(VX, VY, EToV, tol, x0, y0, L1, L2):
    upper_points = set()
    lower_points = set()
    for e, indices in EToV.items():
        for r in range(3):
            s = (r + 1) % 3
            v_r = indices[r]
            v_s = indices[s]

            x_r, y_r = VX[v_r - 1], VY[v_r - 1]
            x_s, y_s = VX[v_s - 1], VY[v_s - 1]

            if (np.allclose(x_r, x0, atol=tol) or np.allclose(y_r, y0 + L2, atol=tol)) \
                    and (np.allclose(x_s, x0, atol=tol) or np.allclose(y_s, y0 + L2, atol=tol)):
                lower_points.add((e, r))

            if (np.allclose(y_r, y0, atol=tol) or np.allclose(x_r, x0 + L1, atol=tol)) \
                    and (np.allclose(y_s, y0, atol=tol) or np.allclose(x_s, x0 + L1, atol=tol)):
                upper_points.add((e, r))

    # for idx, overlap in enumerate(upper_points.intersection(lower_points)):
    #     if idx % 2 == 0:
    #         upper_points.remove(overlap)
    #     else:
    #         lower_points.remove(overlap)

    print()

    return list(lower_points), list(upper_points)


def q1q2(x_i, x_j, y_i, y_j, lam1, lam2, test_case, n, k, VX, VY, etov,exercise):
    if exercise=="2_7":
        if test_case == 1:  # q(x,y) is constant
            u_x = 3  # derivative of (3x + 5y - 7) over x
            u_y = 5  # derivative of (3x + 5y - 7) over y
            n = outernormal(n, k + 1, VX, VY, etov)
            n1, n2 = n[0], n[1]
            q = -lam1 * u_x * n1 - lam2 * u_y * n2
            q_1 = q_2 = q / 2 * np.sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2)
            return q_1, q_2
        if test_case == 2:  # q(x,y) is NOT constant
            x_c = np.mean([x_i, x_j])  # (2.42)
            y_c = np.mean([y_i, y_j])
            n = outernormal(n, k + 1, VX, VY, etov)
            n1, n2 = n[0], n[1]
            u_x = np.cos(x_c) * np.sin(y_c)  # derivative of sin(x)*sin(y) over x
            u_y = np.cos(y_c) * np.sin(x_c)  # derivative of sin(x)*sin(y) over y
            q = -lam1 * u_x * n1 - lam2 * u_y * n2
            q_1 = q_2 = q / 2 * np.sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2)
            return q_1, q_2
    elif exercise=="2_8":
        x_c = np.mean([x_i, x_j])  # (2.42)
        y_c = np.mean([y_i, y_j])
        n = outernormal(n, k + 1, VX, VY, etov)
        n1, n2 = n[0], n[1]
        u_x = -np.sin(np.pi *x_c) * np.cos(np.pi*y_c)*np.pi  # derivative cos(πx) cos(πy)  over x
        u_y =- np.cos(np.pi*x_c) * np.sin(np.pi*y_c)*np.pi  # derivative of cos(πx) cos(πy) over y
        q = -lam1 * u_x * n1 - lam2 * u_y * n2
        q_1 = q_2 = q / 2 * np.sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2)
        return q_1, q_2
    raise ValueError("Unknown test case")


def neumann_boundary_conditions(VX, VY, lam1, lam2, EToV, boundary_edges, qt, b1, test_case,exercise):
    b = np.array(b1)
    for e, r in boundary_edges:
        s = (r + 1) % 3
        vertice_tuple = EToV[e]
        global_r = vertice_tuple[r]
        global_s = vertice_tuple[s]
        x_r, y_r = VX[global_r - 1], VY[global_r - 1]
        x_s, y_s = VX[global_s - 1], VY[global_s - 1]
        q1, q2 = q1q2(x_r, x_s, y_r, y_s, lam1, lam2, test_case, e, r, VX, VY, EToV,exercise=exercise)
        b[global_r - 1] += q1  # TODO: validate
        b[global_s - 1] += q2
    return b




def construct_qt_2_7(etov, VX, VY, test_case, lam1=1, lam2=1):
    qt_dict = dict()
    for e, (v1, v2, v3) in etov.items():
        if test_case == 1:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            def f():
                u_xx = 0
                u_yy = 0
                return - u_xx - u_yy

            q1 = f()
            q2 = f()
            q3 = f()

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
        elif test_case == 2:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            def f(x, y):
                u_xx = -np.sin(x) * np.sin(y)
                u_yy = -np.sin(x) * np.sin(y)
                return -u_xx - u_yy

            q1 = f(x1, y1)
            q2 = f(x2, y2)
            q3 = f(x3, y3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
    return qt_dict


def test_case_data_ex_2_7(test_case):
    if test_case == 1:
        x0 = -2.5
        y0 = -4.8
        L1 = 7.6
        L2 = 5.9
        noelms1 = 4
        noelms2 = 3
        X, Y = xy(x0=x0, y0=y0, L1=L1, L2=L2, noelms1=noelms1, noelms2=noelms2)
        etov, M = construct_element_table(noelms1, noelms2)
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
    # nr_of_test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2 = test_case_data_ex_2_7(1)
    nr_of_test_case, X, Y, L1, L2, x0, y0, etov, M, lam1, lam2, qt = test_case_data(2)
    lower_points, upper_points = construct_boundary_edges(X, Y, etov, tol=0.0005, x0=x0, y0=y0, L1=L1, L2=L2)
    visualise(X, Y, etov, lower_points)
    visualise(X, Y, etov, upper_points)
    A, b = assembly(X, Y, etov, lam1=lam1, lam2=lam2, qt=qt, M=M)
    b_neumann1 = neumann_boundary_conditions(X, Y, lam1, lam2, etov, upper_points, qt, b, test_case=nr_of_test_case,exercise="2_7")
    # b_neumann2 = neumann_boundary_conditions(X, Y, lam1, lam2, etov, lower_points, qt, b_neumann1,
    #                                          test_case=nr_of_test_case)
    u = np.linalg.solve(A, b_neumann1)
    print()


if __name__ == '__main__':
    main()
