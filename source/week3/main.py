from functools import lru_cache

from source.week3.construct_etov import construct_element_table
from source.week3.visualisation import plot_2d_grid, plot_surface
from source.week3.xy import xy
import numpy as np


def solve_elements_plane(element_points, elementU):
    A_inv = np.linalg.pinv(element_points)
    plane_params = np.dot(A_inv, elementU)
    return plane_params


@lru_cache(maxsize=10000)
def area(x1, y1, x2, y2, x3, y3):
    S = abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)
    return S


def is_point_in_triangle(x1, y1, x2, y2, x3, y3, x, y):
    A = area(x1, y1, x2, y2, x3, y3)
    A1 = area(x, y, x2, y2, x3, y3)
    A2 = area(x1, y1, x, y, x3, y3)
    A3 = area(x1, y1, x2, y2, x, y)
    decision = np.isclose(A, A1 + A2 + A3, rtol=0.001)
    return decision


def find_element(x_disp, y_disp, X, Y, EToV):
    P = (x_disp, y_disp)
    for e, vertices in EToV.items():
        r, s, t = vertices
        r -= 1
        s -= 1
        t -= 1
        x_r, x_s, x_t = X[r], X[s], X[t]
        y_r, y_s, y_t = Y[r], Y[s], Y[t]
        R = (x_r, y_r)
        S = (x_s, y_s)
        T = (x_t, y_t)

        if is_point_in_triangle(*R, *S, *T, *P):
            return e
    plot_2d_grid(X, Y, EToV, [(x_disp, y_disp)], elements_to_plot=[x for x in EToV.keys()])
    print()


def interpolate_in_mesh(x_val, y_val, element_idx, EToV, X, Y, U):
    r, s, t = EToV[element_idx]
    r -= 1
    s -= 1
    t -= 1
    x_r, x_s, x_t = X[r], X[s], X[t]
    y_r, y_s, y_t = Y[r], Y[s], Y[t]
    elements_matrix = np.array([[x_r, y_r, 1],
                                [x_s, y_s, 1],
                                [x_t, y_t, 1]])
    U_vector = np.array([U[r], U[s], U[t]])
    plane = solve_elements_plane(elements_matrix, U_vector)
    u_hat = np.array([x_val, y_val, 1]).dot(plane)
    return u_hat


def interpolate_in_2d(EToV, U, X, Y, X_disp, Y_disp):
    """
    return 2D vector of the solution interpolated on mesh used for display
    """
    U_hat = np.zeros_like(X_disp)
    for idx, (x_val, y_val) in enumerate(zip(X_disp, Y_disp)):
        element_index = find_element(x_val, y_val, X, Y, EToV)
        u_val = interpolate_in_mesh(x_val, y_val, element_index, EToV, X, Y, U)
        U_hat[idx] = u_val
    return U_hat


def compute_triangle_error(triangle, plane_params1):
    (x_c, y_c, u_c), (x_1, y_1, u_1), (x_2, y_2, u_2) = triangle
    # Find plane P1
    elements_matrix = np.array([[x_1, y_1, 1],
                                [x_2, y_2, 1],
                                [x_c, y_c, 1]])
    U_vector = np.array([u_1, u_2, u_c])

    plane_params2 = solve_elements_plane(elements_matrix, U_vector)
    # print("Now we have all we need to compute the integrals")
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


def compute_error(element_idx, EToV, X, Y, U_function):
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

    U_r = U_function(x_r, y_r)
    U_s = U_function(x_s, y_s)
    U_t = U_function(x_t, y_t)
    U_c = U_function(x_c, y_c)  # This should be replaced with FEM once we have it

    R = (x_r, y_r, U_r)
    S = (x_s, y_s, U_s)
    T = (x_t, y_t, U_t)
    C = (x_c, y_c, U_c)

    triangles = [(C, R, S), (C, S, T), (C, T, R)]
    elements_matrix = np.array([[x_r, y_r, 1],
                                [x_s, y_s, 1],
                                [x_t, y_t, 1]])
    U_vector = np.array([U_r, U_s, U_t])
    common_plane_params = solve_elements_plane(elements_matrix, U_vector)

    error_sum = sum(compute_triangle_error(triangle, common_plane_params) for triangle in triangles)
    return error_sum


def U_true(X, Y):
    return np.exp(X) + np.exp(Y)


def main():
    elem1 = 4
    elem2 = 3
    x0 = 0
    L1 = 1
    y0 = 0
    L2 = 1

    noelem_display1 = 25
    noelem_display2 = 25

    EToV, M = construct_element_table(elem1, elem2)

    X, Y = xy(x0, y0, L1, L2, elem1, elem2, as_list=True)
    for i in range(500):
        errors = np.array([compute_error(e + 1, EToV, X, Y, U_true) for e in range(len(X))])
        argmax = np.argmax(errors) + 1
        refine_mesh(argmax, EToV, X, Y, U_true)

    X_display, Y_display = xy(x0, y0, L1, L2, noelem_display1, noelem_display2)
    U = U_true(X, Y)
    U_hat = interpolate_in_2d(EToV, U, X, Y, X_display, Y_display)
    plot_surface(X_display, Y_display, noelem_display1, noelem_display2, U_hat)
    plot_2d_grid(X, Y, EToV)


if __name__ == '__main__':
    main()
