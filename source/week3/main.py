from functools import lru_cache

from source.week2.ex2_2a import basfun
from source.week2.ex2_3 import compute_k_rsn, assembly
from source.week3.construct_etov import construct_element_table
from source.week3.integration import integration
from source.week3.mesh_division import newest_node_bisection
from source.week3.visualisation import plot_2d_grid, plot_triangulation, plot_triangulation_old
from source.week3.xy import xy
import numpy as np


def solve_elements_plane(x0, x1, x2, y0, y1, y2, z0, z1, z2):
    ux, uy, uz = [x1 - x0, y1 - y0, z1 - z0]
    vx, vy, vz = [x2 - x0, y2 - y0, z2 - z0]

    u_cross_v = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]
    point = np.array([x0, y0, z0])
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


c_points = dict()


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


def compute_error(element_idx, EToV, X, Y, U_function, element_to_base):
    """
    Computes L2 error that could be reduced by introducing another point into mesh triangle in its centroid.
    :return: float
    """
    r, s, t = EToV[element_idx]
    pt1, pt2 = element_to_base[element_idx]

    r -= 1
    s -= 1
    t -= 1
    pt1 -= 1
    pt2 -= 1

    x_r, x_s, x_t, x_1, x_2 = X[r], X[s], X[t], X[pt1], X[pt2]
    y_r, y_s, y_t, y_1, y_2 = Y[r], Y[s], Y[t], Y[pt1], Y[pt2]

    x_c = (x_1 + x_2) / 2
    y_c = (y_1 + y_2) / 2

    u_r = U_function(x_r, y_r)
    u_s = U_function(x_s, y_s)
    u_t = U_function(x_t, y_t)

    common_plane = solve_elements_plane(x_r, x_s, x_t, y_r, y_s, y_t, u_r, u_s, u_t)

    u_c = eval_u_on_plane(common_plane, x_c, y_c)  # This should be replaced with FEM once we have it
    u_d = U_true(x_c, y_c)
    diff = np.abs(u_d - u_c)
    return diff


def q_tilda(x, y):
    u_xx = np.exp(-100 * (-0.5 + x) ** 2 + (-0.75 + y) ** 2) * (9800 - 40000 * x + 40000 * x ** 2)
    u_yy = np.exp(-100 * (-0.5 + x) ** 2 + (-0.75 + y) ** 2) * (4.25 - 6 * y + 4 * y ** 2)
    return -u_xx - u_yy  # -q_tilda(x,y)=u_xx + u_yy, where u(x,y)=np.exp(-100 * ((X - 0.5) ** 2 + (Y - 0.75) ** 2))


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
        return np.exp(-100 * ((x - 0.5) ** 2 + (y - 0.75) ** 2))  # f(x,y)=u(x,y) in our test case

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
    return np.exp(-100 * ((X - 0.5) ** 2 + (Y - 0.75) ** 2))


def U_hat(X, Y, EToV):
    L1 = 1
    L2 = 1
    x0 = 0
    y0 = 0
    qt = construct_qt(EToV, X, Y)
    N = len(X)
    A, b = assembly(X, Y, EToV, lam1=1, lam2=1, qt=qt, M=N)
    A_updated, b_updated = boundary_conditions(X, Y, L1=L1, L2=L2, x0=x0, y0=y0, A=A, b=b)
    u_hat = np.linalg.solve(A_updated, b_updated)
    return u_hat


def main():
    elem1 = 3
    elem2 = 3
    L1 = 1
    L2 = 1
    x0 = 0
    y0 = 0

    EToV, element_to_base, base_to_elements, edge_to_element, M = construct_element_table(elem1, elem2)

    X, Y = xy(x0, y0, L1, L2, elem1, elem2, as_list=True)
    plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    optimization_steps = 0
    tol = 0.01
    max_error = tol + 1

    while max_error > tol:
        errors = np.array([compute_error(e, EToV, X, Y, U_true, element_to_base)
                           for e in range(1, len(EToV) + 1)])
        argmax = np.argmax(errors)
        max_error = errors[argmax]
        print(max_error)
        EToV, X, Y, element_to_base, edge_to_element = newest_node_bisection(argmax + 1, EToV, X, Y, element_to_base,
                                                                             edge_to_element)
        optimization_steps += 1

        # plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    plot_2d_grid(X, Y, EToV, text=False)
    plot_triangulation_old(EToV, X, Y, U_true)


if __name__ == '__main__':
    main()
