import math

import numpy as np
from ex2_1 import construct_element_table, xy
from ex2_3 import assembly
from ex2_4 import boundary_conditions
from ex2_5 import u_hat, max_error
from source.week2.ex2_6 import neumann_boundary_conditions, construct_boundary_edges
from source.week2.ex2_7 import dirichlet_boundary_conditions


def u(x, y):
    return math.cos(math.pi * x) * math.cos(math.pi * y)


def neumann_boundary():
    pass


def construct_qt(etov_dict, VX, VY, test_case=None):
    """
    Computes points as means of each elements' vertices.
    """
    qt_dict = dict()
    for e, (v1, v2, v3) in etov_dict.items():
        if test_case == 1:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            def f(x, y):
                return -2 * math.pi ** 2 * math.cos(math.pi * x) * math.cos(math.pi * y)

            q1 = f(x1, y1)
            q2 = f(x2, y2)
            q3 = f(x3, y3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
        elif test_case == 2:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            def f(x, y):
                return -2 * y ** 2 - 2 * x ** 2

            q1 = f(x1, y1)
            q2 = f(x2, y2)
            q3 = f(x3, y3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
        else:
            raise Exception("Unknown test case")
    return qt_dict


def conditions(part=None):
    if part == "b":
        L1 = 1
        L2 = 1
        x0 = 0
        y0 = 0
        lam1 = 1
        lam2 = 1
        noelms1 = 3
        noelms2 = 3
        return L1, L2, x0, y0, lam1, lam2, noelms1, noelms2

    if part == "c":
        L1 = 1
        L2 = 1
        x0 = -1
        y0 = -1
        lam1 = 1
        lam2 = 1
        noelms1 = 6
        noelms2 = 6
        return L1, L2, x0, y0, lam1, lam2, noelms1, noelms2


def main():
    part = "b"
    if part == "b":
        # part b) x0=y0=0,L1=L2=1 --> 0<=x,y<=1
        L1_b, L2_b, x0_b, y0_b, lam1_b, lam2_b, noelms1_b, noelms2_b = conditions(part)
        X_b, Y_b = xy(x0=x0_b, y0=y0_b, L1=L1_b, L2=L2_b, noelms1=noelms1_b, noelms2=noelms2_b)
        etov_dict_b, M_b = construct_element_table(noelms1_b, noelms2_b)
        qt_b = construct_qt(etov_dict_b, X_b, Y_b, test_case=1)
        A_b, b_b = assembly(X_b, Y_b, etov_dict_b, lam1=lam1_b, lam2=lam2_b, qt=qt_b, M=M_b)
        print(np.linalg.det(A_b))
        upper_edges, lower_edges = construct_boundary_edges(X_b, Y_b, etov_dict_b, tol=0.0005, x0=x0_b, y0=y0_b, L1=L1_b, L2=L2_b)
        b_neumann1 = neumann_boundary_conditions(X_b, Y_b, lam1_b, lam2_b, etov_dict_b, upper_edges, qt_b, b_b, test_case=1,exercise="2_8")
        A_final, b_final = dirichlet_boundary_conditions(A_b, b_neumann1, lower_edges, etov_dict_b, X_b, Y_b, "2_8",
                                                         test_case=1)
        print(np.linalg.det(A_final))
        uhat_b = np.linalg.solve(A_final, b_final)

        print("u_hat(x,y) part b): ", uhat_b)
        u_b = []
        for x, y in zip(X_b, Y_b):
            u_b.append(u(x, y))
        print("u(x,y) part b): ", u_b)
        error_b = max_error(u_b, uhat_b, M_b)
        print("Max error part b): ", error_b)

    elif part == "c":
        # part c)
        L1_c, L2_c, x0_c, y0_c, lam1_c, lam2_c, noelms1_c, noelms2_c = conditions(part)
        X_c, Y_c = xy(x0=x0_c, y0=y0_c, L1=L1_c, L2=L2_c, noelms1=noelms1_c, noelms2=noelms2_c)
        etov_dict_c, M_c = construct_element_table(noelms1_c, noelms2_c)
        qt_c = construct_qt(etov_dict_c, X_c, Y_c, test_case=1)
        A_c, b_c = assembly(X_c, Y_c, etov_dict_c, lam1=lam1_c, lam2=lam2_c, qt=qt_c, M=M_c)
        uhat_c = u_hat(X_c, Y_c, L1_c, L2_c, x0_c, y0_c, M_c, A_c, b_c, exercise="2_8")
        print("u_hat(x,y) part b): ", uhat_c)
        u_c = []
        for x, y in zip(X_c, Y_c):
            u_c.append(u(x, y))
        print("u(x,y) part b): ", u_c)
        error_c = max_error(u_c, uhat_c, M_c)
        print("Max error part b): ", error_c)


if __name__ == '__main__':
    main()
