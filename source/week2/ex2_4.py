from collections import defaultdict

import numpy as np
from ex2_1 import construct_element_table, xy
from ex2_3 import assembly, construct_qt

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


def f(x, y):
    return x ** 3 - x ** 2 * y + y ** 2 - 1


def boundary_conditions(X, Y, EToV, L1, L2, x0, y0, M, A, b, test_case=None, b2_target=None):
    # which_values_are_correct = [True, True, True, True, False, False, False, True, False, False, False,
    #                             True, False, False, False, True, False, False, False, True]

    which_values_are_correct = [True, True, True, True, True, False, False, True, True, False, False, True, True, False,
                                False, True, True, True, True, True]
    b_updates = defaultdict(list)

    for idx, val in enumerate(b):
        b_updates[idx].append(val)
    for i in range(len(X)):
        if np.allclose(X[i], x0) or np.allclose(X[i], x0 + L1) or np.allclose(Y[i], y0) or np.allclose(Y[i], y0 + L2):

            A[i, i] = 1

            if test_case == 1:
                is_this_update_correct1 = which_values_are_correct[i]
                b[i] = 1
                b_updates[i].append(f"=1")
            elif test_case == 2:
                is_this_update_correct1 = which_values_are_correct[i]
                b[i] = f(X[i], Y[i])
                b_updates[i].append(f"=f({f(X[i], Y[i])})")

            for j in range(0, M):
                if j != i:
                    A[i, j] = 0
                    if x0 < X[j] < x0 + L1 and y0 < Y[j] < y0 + L2:
                        is_this_update_correct2 = which_values_are_correct[j]

                        b[j] -= A[j, i] * f(X[i], Y[i])
                        b_updates[j].append((f"-={A[j, i]} * {b[i]} [{A[j, i] * b[i]}]", i, j))

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
    ntest_case = 2
    nr_of_test_case, X, Y, L1, L2, x0, y0, etov_dict, M, lam1, lam2, qt = test_case_data(ntest_case)

    b2_target = [-22.2900, -10.4572, 9.1111, 36.4150, -0.4020, -11.7107, 16.9886, 23.5520, 0.5480, -42.5668, -26.5414,
                 32.3490, 21.7140, 95.6645, 220.9047, 103.9600, 104.2500, 154.9441, 213.3738, 279.5390]

    A, b = assembly(X, Y, etov_dict, lam1=lam1, lam2=lam2, qt=qt, M=M)
    b_target = [9.3832, 9.9506, 2.6018, 0.5674, -0.8996, -0.4982, -15.1958, -6.9474, -22.1986, -43.0962, -57.7938,
                -28.2464, -43.4976, -85.6942, -100.3918, -49.5454, -17.3824,
                -53.0468, -60.3956, -43.0132, ]

    A_updated, b_updated = boundary_conditions(X, Y, etov_dict, L1=L1, L2=L2, x0=x0, y0=y0, M=M, A=A, b=b,
                                               test_case=ntest_case, b2_target=b2_target)
    # b2_target = [1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 2.0833, 2.0833, 1.0000, 1.0000, 0.7500, 0.7500, 1.0000, 1.0000,
    #              2.0833, 2.0833, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000]
    updated_b_indices = b == b_updated
    which_b_values_are_correct = np.isclose(b_updated, b2_target)
    print()


if __name__ == '__main__':
    main()
