from functools import lru_cache

from source.week2.ex2_2a import basfun
from source.week2.ex2_3 import compute_k_rsn
from source.week3.construct_etov import construct_element_table
from source.week3.integration import integration
from source.week3.visualisation import plot_2d_grid, plot_surface, plot_triangulation
from source.week3.xy import xy
import numpy as np


@lru_cache(maxsize=10000)
def area(x1, y1, x2, y2, x3, y3):
    S = abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)
    return S


def is_point_in_triangle(x1, y1, x2, y2, x3, y3, x, y):
    A = area(x1, y1, x2, y2, x3, y3)
    A1 = area(x, y, x2, y2, x3, y3)
    A2 = area(x1, y1, x, y, x3, y3)
    A3 = area(x1, y1, x2, y2, x, y)
    decision = np.isclose(A, A1 + A2 + A3, rtol=0.0000001)
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


def interpolate_in_mesh(x_val, y_val, element_idx, EToV, X, Y, U):
    r, s, t = EToV[element_idx]
    r -= 1
    s -= 1
    t -= 1
    x_r, x_s, x_t = X[r], X[s], X[t]
    y_r, y_s, y_t = Y[r], Y[s], Y[t]
    plane = solve_elements_plane(x_r, x_s, x_t, y_r, y_s, y_t, U[r], U[s], U[t])
    u_hat = eval_u_on_plane(plane, x_val, y_val)
    # u = U_true(x_val, y_val)
    # diff = np.abs(u_hat-u)
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

    u_r = U_function(x_r, y_r)
    u_s = U_function(x_s, y_s)
    u_t = U_function(x_t, y_t)
    u_c = U_function(x_c, y_c)  # This should be replaced with FEM once we have it

    common_plane = solve_elements_plane(x_r, x_s, x_t, y_r, y_s, y_t, u_r, u_s, u_t)
    u_d = eval_u_on_plane(common_plane, x_c, y_c)
    diff = np.abs(u_d - u_c)
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


def assembly(VX, VY, EToV, lam1, lam2, qt, M):
    N = len(EToV)
    A = np.zeros((M, M))
    b = np.zeros(M)
    for n in range(1, N + 1):
        abc, delta = basfun(n, VX, VY, EToV)
        for r in range(3):
            i = EToV[n][r] - 1
            qr = qt[n] * np.abs(delta) / 3  # (2.28)
            b[i] += qr
            for s in range(3):
                j = EToV[n][s] - 1
                k_rsn = compute_k_rsn(abc, delta, r, s, lam1, lam2)  # (2.27)
                A[i, j] += k_rsn
    return A, b


def boundary_conditions(X, Y, L1, L2, x0, y0, M, A, b):
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
    elem1 = 5
    elem2 = 4
    x0 = 2.5
    L1 = 7.6
    y0 = 4.8
    L2 = 5.9

    noelem_display1 = 35
    noelem_display2 = 35

    EToV, M = construct_element_table(elem1, elem2)

    X, Y = xy(x0, y0, L1, L2, elem1, elem2, as_list=True)
    plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    optimization_steps = 0
    tol = 0.01
    max_error = tol + 1
    qt = construct_qt(EToV, X, Y)
    A, b = assembly(X, Y, EToV, lam1=1, lam2=1, qt=qt, M=M)
    A_updated, b_updated = boundary_conditions(X, Y, L1=L1, L2=L2, x0=x0, y0=y0, M=M, A=A, b=b)
    u_hat = np.linalg.solve(A_updated, b_updated)
    u=U_true(X,Y)
    print("u_hat",u_hat)
    print("u", u)

    while max_error > tol:
        qt = construct_qt(EToV, X, Y)
        A, b = assembly(X, Y, EToV, lam1=1, lam2=1, qt=qt, M=M)
        A_updated, b_updated = boundary_conditions(X, Y, L1=L1, L2=L2, x0=x0, y0=y0, M=M, A=A, b=b)
        u_hat = np.linalg.solve(A_updated, b_updated)
        # plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))
        errors = np.array([compute_error(e, EToV, X, Y, U_true) for e in range(1, len(EToV) + 1)])
        argmax = np.argmax(errors) + 1
        max_error = errors[argmax]
        EToV, X, Y = refine_mesh(argmax, EToV, X, Y, U_true)

        optimization_steps += 1
    print(optimization_steps, max_error)
    plot_2d_grid(X, Y, EToV, text=False)
    plot_triangulation(EToV, X, Y, U_true)


if __name__ == '__main__':
    main()
