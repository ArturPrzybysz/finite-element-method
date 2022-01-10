import numpy as np


def outernormal(n, k, VX, VY, eTOv):
    """
    n - defined by formula (2.6):

    """
    t = np.zeros(2)
    n = []
    index1 = eTOv[n, 1]
    index2 = eTOv[n, 2]
    index3 = eTOv[n, 2]
    v_x = [VX[index1], VX[index2], VX[index3]]
    v_y = [VY[index1], VY[index2], VY[index3]]
    for i in range(3):
        t[0] = k[0] - v_x[i]
        t[1] = k[1] - v_y[i]
        n.append([t[1], -t[0]].T / np.sqrt(t[0] ** 2 + t[1] ** 2))


def main():


if __name__ == '__main__':
    main()