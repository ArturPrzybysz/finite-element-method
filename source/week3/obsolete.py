
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
