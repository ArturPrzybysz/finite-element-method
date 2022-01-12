import numpy as np
from source.week2.ex2_1 import construct_element_table, xy
from source.week2.ex2_3 import assembly, construct_qt
from source.week2.ex2_4 import boundary_conditions



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

    A, b = assembly(X, Y, etov_dict, lam1=lam1, lam2=lam2, qt=qt, M=M)

    A, b = boundary_conditions(X, Y, etov_dict, L1=L1, L2=L2, x0=x0, y0=y0, M=M, A=A, b=b, test_case=ntest_case)


if __name__ == '__main__':
    main()