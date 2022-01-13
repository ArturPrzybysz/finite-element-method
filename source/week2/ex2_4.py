import numpy as np
from ex2_1 import construct_element_table, xy
from ex2_3 import assembly, construct_qt, sparse_diags

'''
def somega(EToV, element_idx, M):
    i_idx, j, k = EToV[element_idx]
    somega = []
    omega = []
    for i in range(1, (M+1)):
        i_tocheck, j_tocheck, k_tocheck = EToV[i]
        if (i_idx == j_tocheck):
            omega = omega.append(i_tocheck, k_tocheck)
        if (i_idx == k_tocheck):
           omega = omega.append(i_tocheck, j_tocheck)
        somega = somega.append(omega)   
    return somega
'''

"""
def f(x, y, test_case=None):
    if test_case == 1:
        return 1
    if test_case == 2:
        return x ** 3 - x ** 2 * y + y ** 2 - 1



def f(x, y, test_case=None):
    if test_case == 1:
        return x ** 3 - x ** 2 * y + y ** 2 - 1
    if test_case == 2:
        return x ** 2 * y ** 2

"""


def boundary_conditions(X, Y, L1, L2, x0, y0, M, A, b, test_case=None, exercise=None):
    if exercise == "2_4":
        def f(x, y, test_case=None):
            if test_case == 1:
                return 1
            if test_case == 2:
                return x ** 3 - x ** 2 * y + y ** 2 - 1
    elif exercise == "2_5":
        def f(x, y, test_case=None):
            if test_case == 1:
                return x ** 3 - x ** 2 * y + y ** 2 - 1
            if test_case == 2:
                return x ** 2 * y ** 2
    else:
        raise Exception("Unknown exercise")

    for i in range(len(X)):
        if np.allclose(X[i], x0) or np.allclose(X[i], x0 + L1) or np.allclose(Y[i], y0) or np.allclose(Y[i], y0 + L2):
            A[i, i] = 1
            b[i] = f(X[i], Y[i], test_case)

            for j in range(0, M):
                if j != i:
                    A[i, j] = 0
                    if x0 < X[j] < x0 + L1 and y0 < Y[j] < y0 + L2:
                        b[j] -= A[j, i] * f(X[i], Y[i], test_case)
                        A[j, i] = 0

    return np.round(A, 4), np.round(b, 4)


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


def main():
    ntest_case = 1
    nr_of_test_case, X, Y, L1, L2, x0, y0, etov_dict, M, lam1, lam2, qt = test_case_data(ntest_case)

    A, b = assembly(X, Y, etov_dict, lam1=lam1, lam2=lam2, qt=qt, M=M)
    A_updated, b_updated = boundary_conditions(X, Y, L1=L1, L2=L2, x0=x0, y0=y0, M=M, A=A, b=b,
                                               test_case=ntest_case, exercise="2_4")
    B, d = sparse_diags(A_updated)
    print("d", d)
    print("b: ", b)


if __name__ == '__main__':
    main()
