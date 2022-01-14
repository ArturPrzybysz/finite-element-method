import numpy as np

from source.week2.ex2_1 import construct_element_table, xy


def outernormal(n, k, VX, VY, EtoV):
    """
    n - index of an element in EtoV dict
    k - number of an edge: only in {1, 2, 3}
    VX - array of X points values
    VY - array of Y points values
    EtoV - dictionary connecting element index with node indices
    """

    v1, v2, v3 = EtoV[n]

    v_x = np.array([VX[v1 - 1], VX[v2 - 1], VX[v3 - 1]])
    v_y = np.array([VY[v1 - 1], VY[v2 - 1], VY[v3 - 1]])

    point1, point2 = None, None
    if k == 1:
        point1 = np.array([v_x[0], v_y[0]])
        point2 = np.array([v_x[1], v_y[1]])
    if k == 2:
        point1 = np.array([v_x[1], v_y[1]])
        point2 = np.array([v_x[2], v_y[2]])
    if k == 3:
        point1 = np.array([v_x[2], v_y[2]])
        point2 = np.array([v_x[0], v_y[0]])

    # Formula (2.6)
    t12 = - np.array([point2[0] - point1[0], point2[1] - point1[1]])
    numerator = np.array([t12[1], -t12[0]]).T
    denominator = np.sqrt(t12[1] ** 2 + t12[0] ** 2)
    n = numerator / denominator

    # Validation
    test = n.dot(t12)
    assert test == 0
    return n


def main():
    EtoV = construct_element_table(4, 3)
    VX, VY = xy(-2.5, -4.8, 7.6, 5.9, 4, 3)

    outernormal(9, k=3, VX=VX, VY=VY, EtoV=EtoV)


if __name__ == '__main__':
    main()
