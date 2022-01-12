"""
a) Write a Matlab routine that constructs an array for holding information about boundary
edges (beds) of the triangular mesh. (HINT: the use of a distance function fd can make the
job of finding nodes that belong to a boundary implicitly defined by the distance function
very easy. Below it is included in the argument list as an optional argument. top is a
tolerance to be used with the signed distance function to tune the selection of nodes.)
"""
import numpy as np

from source.week2.ex2_3 import test_case_data

"""
Hint: Check slide 30 in lecture 4 on Mesh Generation. Signed distance functions are available in the DistMesh package
"""


def construct_beds(VX, VY, EToV, tol, u_function):
    """

    """
    Ms = []
    rs = []
    for e, indices in EToV.items():
        for r, v in enumerate(indices):
            x, y = VX[v - 1], VY[v - 1]
            y_u = u_function(x)
            if np.linalg.norm(np.array([x, y]) - np.array(x, y_u)) < tol:
                Ms.append(e)
                rs.append(r)
    return np.array([Ms, rs])


def u_func1(x):


def main():
    nr_of_test_case, X, Y, etov, M, lam1, lam2, qt = test_case_data(1)
    construct_beds(X, Y, etov, tol=0.000005, u_func1)


if __name__ == '__main__':
    main()
