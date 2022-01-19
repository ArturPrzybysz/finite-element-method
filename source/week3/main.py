from functools import lru_cache

from source.week2.ex2_2a import basfun
from source.week2.ex2_3 import compute_k_rsn, assembly
from source.week3.construct_etov import construct_element_table
from source.week3.integration import integration
from source.week3.visualisation import plot_2d_grid, plot_surface, plot_triangulation
from source.week3.xy import xy
import numpy as np


def solve_elements_plane(x0, x1, x2, y0, y1, y2, z0, z1, z2):
    ux, uy, uz = [x1 - x0, y1 - y0, z1 - z0]
    vx, vy, vz = [x2 - x0, y2 - y0, z2 - z0]

    u_cross_v = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]
    p0 = x0, y0, z0
    point = np.array(p0)
    normal = np.array(u_cross_v)

    a = normal[0]
    b = normal[1]
    c = normal[2]
    d = -point.dot(normal)
    return a, b, c, d


def eval_u_on_plane(plane, x, y):
    a, b, c, d = plane
    u = (-a * x - b * y - d) / c
    return u


def compute_triangle_error(triangle, common_plane, D):
    (x_c, y_c, u_c), (x_1, y_1, u_1), (x_2, y_2, u_2) = triangle
    # Find plane P1
    plane = solve_elements_plane(x_1, x_2, x_c, y_1, y_2, y_c, u_1, u_2, u_c)
    A1, A2, A3, A4 = common_plane  # lower plane: U1(x, y) = (-A1*x - A2*y - A4) / A3
    A5, A6, A7, A8 = plane  # upper plane: U2(x, y) = (-A5*x - A6*y - A8) / A7

    integral1 = integration(common_plane, plane, )
    integral2 = integration(common_plane, plane, )
    result = integral1 + integral2
    return np.random.rand()


def refine_mesh(element_idx, EToV, X, Y, U_function):
    r, s, t = EToV[element_idx]
    r -= 1
    s -= 1
    t -= 1
    x_r, x_s, x_t = X[r], X[s], X[t]
    y_r, y_s, y_t = Y[r], Y[s], Y[t]
    x_c = (x_r + x_s + x_t) / 3  # Point that splits currently considered mesh into 3 parts (smaller triangles)
    y_c = (y_r + y_s + y_t) / 3
    X.append(x_c)
    Y.append(y_c)
    c = len(X)
    triangle1 = (r + 1, s + 1, c)
    triangle2 = (s + 1, t + 1, c)
    triangle3 = (t + 1, r + 1, c)
    EToV[element_idx] = triangle1
    last_idx = len(EToV)
    EToV[last_idx + 1] = triangle2
    EToV[last_idx + 2] = triangle3
    return EToV, X, Y


def compute_error(element_idx, EToV, X, Y, U_function, last_error):
    """
    Computes L2 error that could be reduced by introducing another point into mesh triangle in its centroid.
    :return: float
    """
    r, s, t = EToV[element_idx]
    r -= 1
    s -= 1
    t -= 1
    x_r, x_s, x_t = X[r], X[s], X[t]
    y_r, y_s, y_t = Y[r], Y[s], Y[t]

    x_c = (x_r + x_s + x_t) / 3  # Point that splits currently analysed mesh into 3 parts (smaller triangles)
    y_c = (y_r + y_s + y_t) / 3

    u_r = U_function(x_r, y_r)
    u_s = U_function(x_s, y_s)
    u_t = U_function(x_t, y_t)
    common_plane = solve_elements_plane(x_r, x_s, x_t, y_r, y_s, y_t, u_r, u_s, u_t)

    u_c = eval_u_on_plane(common_plane, x_c, y_c)  # This should be replaced with FEM once we have it

    # u_d = eval_u_on_plane(common_plane, x_c, y_c)
    u_d = U_true(x_c, y_c)
    diff = np.abs(u_d - u_c)
    # for i in range(3):
    #     plane = ...
    #     integration(common_plane,plane, )

    return diff


def q_tilda(x, y):
    return -6 * x + 2 * y - 2  # -q_tilda(x,y)=u_xx + u_yy, where u(x,y)=x^3 -x^2 y +y^2-1


def construct_qt(etov_dict, VX, VY):
    """
    Computes points as means of each elements' vertices.
    """
    qt_dict = dict()
    for e, (v1, v2, v3) in etov_dict.items():
        x1, x2, x3 = VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]
        y1, y2, y3 = VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]

        q1 = q_tilda(x1, y1)
        q2 = q_tilda(x2, y2)
        q3 = q_tilda(x3, y3)

        q_avg = np.mean([q1, q2, q3])
        qt_dict[e] = q_avg

    return qt_dict


def boundary_conditions(X, Y, L1, L2, x0, y0, A, b):
    def f(x, y):
        return x ** 3 - x ** 2 * y + y ** 2 - 1  # f(x,y)=u(x,y) in our test case

    for i in range(len(X)):
        if np.allclose(X[i], x0) or np.allclose(X[i], x0 + L1) or np.allclose(Y[i], y0) or np.allclose(Y[i], y0 + L2):
            A[i, i] = 1
            b[i] = f(X[i], Y[i])

            for j in range(0, len(X)):
                if j != i:
                    A[i, j] = 0
                    if x0 < X[j] < x0 + L1 and y0 < Y[j] < y0 + L2:
                        b[j] -= A[j, i] * f(X[i], Y[i])
                        A[j, i] = 0

    return A, b


def U_true(X, Y):
    X = np.array(X)
    Y = np.array(Y)
    # np.exp(-100 * ((X - 0.5) ** 2 + (Y - 0.75) ** 2))
    return X ** 3 - X ** 2 * Y + Y ** 2 - 1


def main():
    # elem1 = 7
    # elem2 = 7
    # x0 = 0
    # L1 = 1
    # y0 = 0
    # L2 = 1
    elem1 = 4
    elem2 = 3
    L1 = 7.6
    L2 = 5.9
    x0 = -2.5
    y0 = -4.8

    EToV, M = construct_element_table(elem1, elem2)
    X, Y = xy(x0, y0, L1, L2, elem1, elem2, as_list=True)
    # plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    optimization_steps = 0
    tol = 0.0001
    max_error = tol + 1
    last_error = max_error
    qt = construct_qt(EToV, X, Y)
    N = len(X)
    A, b = assembly(X, Y, EToV, lam1=1, lam2=1, qt=qt, M=N)
    A_updated, b_updated = boundary_conditions(X, Y, L1=L1, L2=L2, x0=x0, y0=y0, A=A, b=b)
    u_hat = np.linalg.solve(A_updated, b_updated)
    u = U_true(X, Y)
    print("u_hat", u_hat)
    print("u", u)
    print(max(u_hat - u))


""""
    while max_error > tol:
        errors = np.array([compute_error(e, EToV, X, Y, U_true, last_error) for e in range(1, len(EToV) + 1)])
        argmax = np.argmax(errors)
        qt = construct_qt(EToV, X, Y)
        A, b = assembly(X, Y, EToV, lam1=1, lam2=1, qt=qt, M=M)

        A_updated, b_updated = boundary_conditions(X, Y, L1=L1, L2=L2, x0=x0, y0=y0, M=M, A=A, b=b)
        u_hat = np.linalg.solve(A_updated, b_updated)
        # plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))
        errors = np.array([compute_error(e, EToV, X, Y, U_true) for e in range(1, len(EToV) + 1)])
        argmax = np.argmax(errors) + 1
        max_error = errors[argmax]
        last_error = max_error
        print(max_error)
        EToV, X, Y = refine_mesh(argmax + 1, EToV, X, Y, U_true)
        # plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

        EToV, X, Y = refine_mesh(argmax, EToV, X, Y, U_true)

        optimization_steps += 1
    plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    print(optimization_steps, max_error)
    plot_2d_grid(X, Y, EToV, text=False)
    plot_triangulation(EToV, X, Y, U_true)
"""

if __name__ == '__main__':
    main()
