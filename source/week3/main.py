from functools import lru_cache

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


def U_true(X, Y):
    X = np.array(X)
    Y = np.array(Y)
    return np.exp(-100 * ((X - 0.5) ** 2 + (Y - 0.75) ** 2))


def main():
    elem1 = 7
    elem2 = 7
    x0 = 0
    L1 = 1
    y0 = 0
    L2 = 1

    EToV, M = construct_element_table(elem1, elem2)

    X, Y = xy(x0, y0, L1, L2, elem1, elem2, as_list=True)
    plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    optimization_steps = 0
    tol = 0.0001
    max_error = tol + 1
    last_error =max_error

    while max_error > tol:
        errors = np.array([compute_error(e, EToV, X, Y, U_true, last_error) for e in range(1, len(EToV) + 1)])
        argmax = np.argmax(errors)
        max_error = errors[argmax]
        last_error = max_error
        print(max_error)
        EToV, X, Y = refine_mesh(argmax + 1, EToV, X, Y, U_true)
        # plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

        optimization_steps += 1
    plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    print(optimization_steps, max_error)
    plot_2d_grid(X, Y, EToV, text=False)
    plot_triangulation(EToV, X, Y, U_true)


if __name__ == '__main__':
    main()
