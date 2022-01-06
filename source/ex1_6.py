import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from source.ex1_2 import plot_solution, interpolate


def u_function(x):
    return np.exp(-800 * (x - 0.4) ** 2) + 0.25 * np.exp(-40 * (x - 0.8) ** 2)


# TODO: plot second derivative versus mesh density


def point_on_line(p1, p2, p3):
    l2 = np.sum((p1 - p2) ** 2)
    t = max(0, min(1, np.sum((p3 - p1) * (p2 - p1)) / l2))
    projection = p1 + t * (p2 - p1)
    return projection


def compute_L2_error(base_length, height):
    """
    Computes area of a squared
    :param base_length: float, length of a base of a triangle
    :param height: float, height of a triangle 
    :return: area under parabola: y = ax^2 in range [0, base_length]
    """
    if np.isclose(height, 0): return 0
    if np.isclose(base_length, 0): return 0
    a = height ** 2 / base_length ** 2
    return a / 3 * base_length ** 3


def compute_error_decrease(fun, VX, EToV) -> np.array:
    """
    Computes estimate of possible error decrease for each element in mesh.

    :param fun: Function float -> float
    :param VX: dict from point id to its position on x axis.
    :param EToV: dict from element id to a tuple of its boundary points.
    """
    L2_loss = dict()
    for e, (idx1, idx2) in EToV.items():
        x1 = VX[idx1]
        x2 = VX[idx2]
        y1 = fun(x1)
        y2 = fun(x2)
        x_half = (x1 + x2) / 2
        y_half = fun(x_half)

        P1 = np.array([x1, y1])
        P2 = np.array([x2, y2])
        P_half = np.array([x_half, y_half])
        P_projected = point_on_line(P1, P2, P_half)

        triangle_height = np.sqrt(np.sum((P_projected - P_half) ** 2))
        base_len1 = np.sqrt(np.sum((P_projected - P1) ** 2))
        base_len2 = np.sqrt(np.sum((P_projected - P2) ** 2))

        L2_loss1 = compute_L2_error(base_len1, triangle_height)
        L2_loss2 = compute_L2_error(base_len2, triangle_height)
        L2_loss[e] = L2_loss1 + L2_loss2

    return L2_loss


def refine_marked(EtoV, VX, idxMarked):
    # [EToVfine, VXfine]
    max_e_idx = max(EtoV.keys())
    max_v_idx = max(VX.keys())
    for idx in idxMarked:
        x1_idx, x2_idx = EtoV[idx]
        x1 = VX[x1_idx]
        x2 = VX[x2_idx]
        x_mid = (x1 + x2) / 2

        if x_mid in VX.values():
            continue

        EtoV[idx] = (x1_idx, max_v_idx + 1)
        EtoV[max_e_idx + 1] = (max_v_idx + 1, x2_idx)
        VX[max_v_idx + 1] = x_mid

        max_e_idx += 1
        max_v_idx += 1
    return EtoV, VX


def refine_until_converged(VX, EtoV, fun, tolerance):
    L2_loss = compute_error_decrease(fun, VX, EtoV)
    step = 0
    while max(L2_loss.values()) >= tolerance:
        marked_idx = [k for k, v in sorted(L2_loss.items(), key=lambda t: t[1])[-1:]]
        EToV, VX = refine_marked(EtoV, VX, marked_idx)
        L2_loss = compute_error_decrease(u_function, VX, EtoV)
        step += 1
    return VX, EtoV


def main():
    """
    QUESTION: "How can the metric be computed? Detail this."
    ANSWER: Having three points: A, B and C, which are left, right bounds and the true point in the middle, we can
    estimate the L2 error using an analytical method. Using basic geometry I find the projection of C onto line
    connecting A and B. This allows to find two perpendicular triangles, which are later squared and the area under
    parabola in the range [0, triangle_base] is computed.
    """
    """
    NOTE1: Even if tolerance is set very low, the function could omit some areas! For example if initialized with 
        3 points, then the left local (and global at the same time) maximum will be omitted. Only when more points
        are used for initialization, this area is 'explored' more than the right local maximum.  
        Show it on examples.
    """
    mesh = np.linspace(start=0, stop=1, num=2)
    u_hat = u_function(mesh)

    x = np.linspace(start=0, stop=1, num=500)
    u_x = u_function(x)
    u_hat = interpolate(u_hat, mesh, x)
    plot_solution(u_hat, u_x, x, mesh)

    VX = {i: x_i for i, x_i in enumerate(mesh)}
    EtoV = {n: (n, n + 1) for n in range(len(mesh) - 1)}

    VX, EtoV = refine_until_converged(VX, EtoV, u_function, tolerance=0.00001)
    mesh2 = np.array(list(sorted([v for v in VX.values()])))
    u_hat = interpolate(u_function(mesh2), mesh2, x)
    plot_solution(u_hat, u_x, x, mesh2)


if __name__ == '__main__':
    main()
