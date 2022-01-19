import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def plot_triangulation(EToV, X, Y, U_true):
    triangles = np.array([v for v in EToV.values()]) - 1
    # triangulation = Triangulation(X, Y, triangles)
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    Z = U_true(X, Y)
    ax.plot_trisurf(X, Y, Z, triangles=triangles)
    plt.show()
    print()

def plot_surface(X, Y, elem1, elem2, U_hat):
    Z_2D = U_hat.reshape((elem2 + 1, elem1 + 1))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = X[::elem2 + 1]
    Y = Y[:elem2 + 1]
    X, Y = np.meshgrid(X, Y)
    surf = ax.plot_surface(X, Y, Z_2D,
                           cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    #
    # # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


def plot_2d_grid(VX, VY, EToV, points_to_plot=[], elements_to_plot=[], heatmap=None, text=True):
    # if heatmap is not None:
    # U_true, X_display, Y_display, elem1, elem2 = heatmap
    # U = U_true(X_display, Y_display)
    # X = X_display[::elem1 + 1]
    # Y = Y_display[:elem1 + 1]
    # X, Y = np.meshgrid(X, Y)
    # A = U.reshape((elem2 + 1, elem1 + 1))
    # plt.pcolormesh(Y, X, A)

    # plt.imshow(A, cmap='hot', interpolation='nearest')
    plt.gca().set_aspect('equal', adjustable='box')

    plt.scatter(VX, VY, c='grey')
    if text:
        for idx, (x, y) in enumerate(zip(VX, VY)):
            plt.text(x, y + 0.01, idx, c="grey")
    for p in points_to_plot:
        plt.scatter(*p, c='red')
    for e in elements_to_plot:
        for r in range(3):
            s = (r + 1) % 3

            v1 = EToV[e][r]
            v2 = EToV[e][s]

            x_r, y_r = VX[v1 - 1], VY[v1 - 1]
            x_s, y_s = VX[v2 - 1], VY[v2 - 1]
            plt.plot([x_r, x_s], [y_r, y_s], c='green', alpha=0.55)
        r, s, t = EToV[e]
        r -= 1
        s -= 1
        t -= 1
        x_r, x_s, x_t = VX[r], VX[s], VX[t]
        y_r, y_s, y_t = VY[r], VY[s], VY[t]
        x_c = (x_r + x_s + x_t) / 3  # Point that splits currently considered mesh into 3 parts (smaller triangles)
        y_c = (y_r + y_s + y_t) / 3
        if text:
            plt.text(x_c, y_c, str(e))

    plt.show()
