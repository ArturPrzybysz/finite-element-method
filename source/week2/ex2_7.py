import numpy as np

from source.week2.ex2_3 import assembly
from source.week2.ex2_6 import test_case_data_ex_2_7, construct_boundary_edges, neumann_boundary_conditions


def u_true(x, y, test_case):
    if test_case == 1:
        return 3 * x + 5 * y - 7
    elif test_case == 2:
        return np.sin(x) * np.sin(y)
    raise Exception()


def main():
    test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2 = test_case_data_ex_2_7(2)
    lower_points, upper_points = construct_boundary_edges(X, Y, etov, tol=0.00005, x0=x0, y0=y0, L1=L1, L2=L2)
    A, b = assembly(X, Y, etov, lam1=lam1, lam2=lam2, qt=qt, M=M)
    b_neumann = neumann_boundary_conditions(X, Y, etov, upper_points, qt, b)
    u_hat = np.linalg.solve(A, b_neumann)
    u = u_true(X, Y, test_case)
    print()


if __name__ == '__main__':
    main()
