import numpy as np


def construct_boundary_edges(VX, VY, EToV, tol, x0, y0, L1, L2):
    upper_points = set()
    lower_points = set()
    lower_nodes, upper_nodes = set(), set()
    for e, indices in EToV.items():
        for r in range(3):
            s = (r + 1) % 3
            v_r = indices[r]
            v_s = indices[s]

            x_r, y_r = VX[v_r - 1], VY[v_r - 1]
            x_s, y_s = VX[v_s - 1], VY[v_s - 1]

            if (np.allclose(x_r, x0, atol=tol) or np.allclose(y_r, y0 + L2, atol=tol)) \
                    and (np.allclose(x_s, x0, atol=tol) or np.allclose(y_s, y0 + L2, atol=tol)):
                lower_points.add((e, r))
                lower_nodes.add(v_r)
                lower_nodes.add(v_s)

            if (np.allclose(y_r, y0, atol=tol) or np.allclose(x_r, x0 + L1, atol=tol)) \
                    and (np.allclose(y_s, y0, atol=tol) or np.allclose(x_s, x0 + L1, atol=tol)):
                upper_points.add((e, r))
                upper_nodes.add(v_r)
                upper_nodes.add(v_s)
    return list(lower_points), list(upper_points), list(lower_nodes), list(upper_nodes)
