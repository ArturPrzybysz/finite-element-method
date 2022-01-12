import numpy as np
from source.week2.ex2_1 import construct_element_table, xy
from source.week2.ex2_3 import assembly, construct_qt
from source.week2.ex2_4 import boundary_conditions
import matplotlib as plt


def u_hat(X, Y, etov_dict, L1, L2, x0, y0, M, A, b, test_case):
    # we have f=u, to get q we need to solve q(x,y)=-(u_(xx)+u_(yy))
    # u_hat --> Au=b

    A, b = boundary_conditions(X, Y, etov_dict, L1=L1, L2=L2, x0=x0, y0=y0, M=M, A=A, b=b, test_case=test_case)
    print("A:", A)
    uhat = np.linalg.solve(A, b)
    return uhat


def u(x, y, test_case=None):
    if test_case == 1:
        return x ** 3 - x ** 2 * y + y ** 2 - 1
    if test_case == 2:
        return x ** 2 * y ** 2


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
                return -6 * x + 2 * y - 2

            q1 = f(x1, y1)
            q2 = f(x2, y2)
            q3 = f(x3, y3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
        elif test_case == 2:
            x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
            y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

            def f(x, y):
                return -4

            q1 = f(x1, y1)
            q2 = f(x2, y2)
            q3 = f(x3, y3)

            q_avg = np.mean([q1, q2, q3])
            qt_dict[e] = q_avg
        else:
            raise Exception("Unknown test case")
    return qt_dict


def test_case_data(nr_of_test_case):
    if nr_of_test_case == 1:
        L1 = 1
        L2 = 1
        x0 = 0
        y0 = 0
        X, Y = xy(x0=x0, y0=y0, L1=L1, L2=L2, noelms1=4, noelms2=3)
        etov_dict, M = construct_element_table(4, 3)
        lam1 = 1
        lam2 = 1
        qt = construct_qt(etov_dict, X, Y, test_case=1)
        return nr_of_test_case, X, Y, L1, L2, x0, y0, etov_dict, M, lam1, lam2, qt

    elif nr_of_test_case == 2:
        L1 = 7.6
        L2 = 5.9
        x0 = -2.5
        y0 = -4.8
        X, Y = xy(x0=x0, y0=y0, L1=L1, L2=L2, noelms1=4, noelms2=3)
        etov_dict, M = construct_element_table(4, 3)
        lam1 = 1
        lam2 = 1
        qt = construct_qt(etov_dict, X, Y, test_case=2)
        return nr_of_test_case, X, Y, L1, L2, x0, y0, etov_dict, M, lam1, lam2, qt
    else:
        raise Exception("Unknown test case")


def max_error(u,u_hat,M):
    errors=[]
    for i in range(M):
        errors.append(abs(u[i]-u_hat[i]))
    return max(errors)


def main():
    ntest_case = 1
    nr_of_test_case, X, Y, L1, L2, x0, y0, etov_dict, M, lam1, lam2, qt = test_case_data(ntest_case)

    A, b = assembly(X, Y, etov_dict, lam1=lam1, lam2=lam2, qt=qt, M=M)
    uhat = u_hat(X, Y, etov_dict, L1, L2, x0, y0, M, A, b, ntest_case)
    print("u_hat", uhat)
    u1 = u(X, Y, ntest_case)
    print("u1", u1)
    error=max_error(u1,uhat,M)
    print(error)

if __name__ == '__main__':
    main()
