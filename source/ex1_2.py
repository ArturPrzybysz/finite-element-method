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
import scipy.linalg
import matplotlib.pyplot as plt


def BVP1D(L: float, c: float, d: float, M: int = None, x: np.array = None, solver="cholesky"):
    """
    Solves 1D bounded value problem.
    Authors: Artur Przybysz, Elena Mongelli, Ivan Knezevic

    Inputs parameters:
    :param L: float, Domain length
    :param c: float, Left boundary condition
    :param d: float, Right boundary condition
    :param M: int: mesh size
    :param x: np.array: 1D mesh vector x(1:{M})
    :param solver: which solver to use ("cholesky", "default")

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
    if solver == "cholesky":
        # L = np.linalg.cholesky(A)
        c, low = scipy.linalg.cho_factor(A)
        u_hat = scipy.linalg.cho_solve((c, low), b)
    else:
        u_hat = np.linalg.solve(A, b)
    return u_hat, x, h


def construct_A_b(h: np.array, c: float, d: float, M: int) -> Tuple[np.array, np.array]:
    """
    Constructs matrix A and vector b using global assembly (Algorithm 1.) and imposes boundary conditions (Algorithm 2.)
    """
    K = np.array([[[1 / h_i + h_i / 3, -1 / h_i + h_i / 6],
                   [-1 / h_i + h_i / 6, [1 / h_i + h_i / 3]]]
                  for h_i in h])

    size = 1 + len(h)
    A = np.zeros((size, size))
    b = np.zeros(size)

    A[0, 0], A[size - 1, size - 1] = 1, 1

    # Iterate over a diagonal of A and fill in values of K (other than ones at 0,0 and L,L)
    for i in range(1, M - 1):
        A[i, i] = K[i - 1, 1, 1] + K[i, 0, 0]
        if i != 1:
            A[i, i - 1] = K[i - 1, 1, 0]
        if i != M - 2:
            A[i, i + 1] = K[i, 0, 1]

    b[0] = c
    b[1] = -K[0, 1, 0] * c
    b[size - 2] = -K[len(K) - 1, 0, 1] * d
    b[size - 1] = d
    return A, b


def plot_solution(u_hat, u, x, mesh):
    plt.plot(x, u, label="$u(x)$")
    plt.plot(x, u_hat, ':', label="$\hat{u}(x)$")
    plt.scatter(mesh, np.zeros_like(mesh))
    plt.legend()
    plt.title('Interpolated vs true $u(x)$.')
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
    diff = np.abs(u - u_hat)
    plt.plot(x, diff)
    plt.title("Difference between $\hat{u}(x)$ and $u(x)$")
    plt.show()


def interpolate(u_hat, mesh, x):
    w = np.zeros_like(x)

    mesh_id = 0
    for i, x_i in enumerate(x):
        if x_i > mesh[mesh_id + 1]:
            mesh_id += 1

        weight1 = x_i - mesh[mesh_id]
        weight2 = mesh[mesh_id + 1] - x_i

        u1 = u_hat[mesh_id]
        u2 = u_hat[mesh_id + 1]
        value = (weight2 * u1 + weight1 * u2) / (weight1 + weight2)
        w[i] = value

    return w


def main():
    L = 2
    c = 1
    d = np.exp(2)
    M = 12

    u_hat, mesh, h = BVP1D(L, c, d, M)
    x = np.linspace(start=min(mesh), stop=max(mesh), num=1000)
    u = u_function(x, a=1, b=0)
    w = interpolate(u_hat, mesh, x)

    plot_solution(w, u, x, mesh)
    plot_difference(w, u, x)
    C = validate(w, u, h)
    print("C=", C)


if __name__ == '__main__':
    main()
