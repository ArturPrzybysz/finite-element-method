"""
Exercise 1.2
The goal of this exercise is to write your first FEM procedure for solving a boundary value
problem in one dimension. This requires that some algorithms of Chapter 1 are implemented.
a) Write a Matlab program based on the skeleton code on the next page that solves the bound-
ary value problem given in Exercise 1.1 using the classical FEM where the local solution is
represented using linear polynomials.

Head:
function [u] = BVP1D(L,c,d,x)

Test case:
Input: L=2, c=1, d=exp(2),
x = [0.0, 0.2, 0.4, 0.6, 0.7, 0.9, 1.4, 1.5, 1.8, 1.9, 2.0]
"""
from typing import Tuple
from time import time_ns
import numpy as np
import sklearn
import matplotlib.pyplot as plt


def BVP1D(L: float, c: float, d: float, M: int = None, x: np.array = None):
    """
    Solves 1D bounded value problem.

    Inputs are following:

    :param L: float, signifies the upper limit of x
    :param c: float, value of u_0
    :param d: float, value of u_{L-1}
    :param M: int: mesh size
    :param x: np.array: mesh given directly to the problem

    :return: u, which is an approximate of solution function
    """

    # generate 1D mesh
    if x is None and M is not None:
        x = np.linspace(0, L, num=M)
    elif x is not None:
        pass
    else:
        raise ValueError("Either M or x must be filled, but neither is given.")

    # lengths of sub-intervals
    h = [x[i + 1] - x[i] for i in range(len(x) - 1)]

    A, b = construct_A_b(h, c, d, M)
    u_hat = np.linalg.solve(A, b)
    return u_hat, x, h


def construct_A_b(h: np.array, c: float, d: float, M: int) -> Tuple[np.array, np.array]:
    """
    Constructs matrix A and vector b which describe a problem
    """
    # Seems wasteful, as K[i, 1, 1] == K[i, 2, 2] and  K[i, 1, 2] == K[i, 2, 1] ???
    K = np.array([[[1 / h_i + h_i / 3, -1 / h_i + h_i / 6],
                   [-1 / h_i + h_i / 6, [1 / h_i + h_i / 3]]]
                  for h_i in h])

    size = 1 + len(h)  # HERE???
    A = np.zeros((size, size))
    b = np.zeros(size)

    A[0, 0], A[size - 1, size - 1] = 1, 1

    # Iterate over a diagonal of A and fill in values of K (other than ones at 0,0 and L,L)
    for i in range(1, M - 1):
        A[i, i] = K[i - 1, 1, 1] + K[i, 0, 0]
        if i != 1:
            A[i, i - 1] = K[i - 1, 1, 0]
        if i != M - 2:  # TODO: fix this
            A[i, i + 1] = K[i, 0, 1]

    b[0] = c
    b[1] = -K[0, 1, 0] * c
    b[size - 2] = -K[len(K) - 1, 0, 1] * d
    b[size - 1] = d
    return A, b


def plot_solution(u_hat, u, x):
    plt.plot(x, u_hat)
    plt.plot(x, u)
    plt.show()


def validate(u_hat, u, h):
    """
    Uses formula (1.33) to compute error.
    """
    C = max(u_hat - u) / max(h) ** 2
    return C


def u_function(x, a=1.0, b=1.0):
    y = a * np.exp(x) + b * np.exp(-x)
    return y


def plot_difference(u_hat, u, x):
    diff = u - u_hat
    plt.plot(x, diff)
    plt.title("Difference between $\hat{u}(x)$ and $u(x)$")
    plt.show()


def main():
    L = 2
    c = 1
    d = np.exp(2)
    M = 12

    u_hat, x, h = BVP1D(L, c, d, M)
    u = u_function(x, a=1, b=0)
    plot_solution(u_hat, u, x)
    plot_difference(u_hat, u, x)
    C = validate(u_hat, u, h)
    print("C=", C)


if __name__ == '__main__':
    main()
