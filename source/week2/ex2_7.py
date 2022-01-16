from itertools import chain
import numpy as np
import matplotlib.pyplot as plt

from source.week2.ex2_3 import assembly
from source.week2.ex2_6 import test_case_data_ex_2_7, construct_boundary_edges, neumann_boundary_conditions, visualise


def u_true(x, y, test_case):
    if test_case == 1:
        return 3 * x + 5 * y - 7
    elif test_case == 2:
        return np.sin(x) * np.sin(y)
    raise Exception()


def find_neighbors(n, etov):
    vertices_with_n = [vertices
                       for e, vertices in etov.items()
                       if n in vertices]
    vertices_set = {v
                    for vertice_tuple in vertices_with_n
                    for v in vertice_tuple}
    vertices_set.remove(n)
    return list(vertices_set)


def find_boundary_indices(lower_edges, etov):
    global_indices_edges = {etov[x[0]][x[1]] for x in lower_edges}
    return list(global_indices_edges)


def f_i_value(exercise, test_case, x, y):
    if exercise == "2_7":
        if test_case == 1:
            return 3 * x + 5 * y - 7
        if test_case == 2:
            return np.sin(x) * np.sin(y)
    raise ValueError("Wrong exercise value.")
    # elif exercise == "2_8":
    #     if test_case == 1:
    #         return x ** 3 - x ** 2 * y + y ** 2 - 1
    #     if test_case == 2:
    #         return x ** 2 * y ** 2



# function [A,b] = dirbc(bnodes,f,A,b)
#     M = length(b);
#     nodes = setdiff(linspace(1,M,M),bnodes);
#     Nodes = zeros(M,1);
#     Nodes(nodes)=1;
#     for i = bnodes
#         A(i,i) = 1;
#         b(i) = f(i);
#         for j = 1:M
#             if i!=j
#                 A(i,j) = 0;
#                 if Nodes(j)
#                     b(j) = b(j)-A(j,i)*f(i);
#                     A(j,i) = 0;
#                 end
#             end
#         end
#     end
# end
def dirichlet_boundary_conditions(A_old, b_old, lower_edge, etov, VX, VY, exercise, test_case):
    A = np.array(A_old)
    b = np.array(b_old)
    S_gamma = find_boundary_indices(lower_edge, etov)

    for i in S_gamma:
        x_i, y_i = VX[i - 1], VY[i - 1]
        f_i = f_i_value(exercise, test_case, x_i, y_i)

        A[i - 1, i - 1] = 1
        b[i - 1] = f_i
        neighbors = find_neighbors(i, etov)

        for j in neighbors:
            if j < i:
                b[j - 1] -= A[j - 1, i - 1] * f_i
                A[j - 1, i - 1] = 0
            else:
                b[j - 1] -= A[i - 1, j - 1] * f_i
                A[i - 1, j - 1] = 0
    return A, b


def main():
    test_case, X, Y, etov, M, lam1, lam2, qt, x0, y0, L1, L2 = test_case_data_ex_2_7(1)
    upper_edges, lower_edges = construct_boundary_edges(X, Y, etov, tol=0.0005, x0=x0, y0=y0, L1=L1, L2=L2)
    visualise(X, Y, etov, upper_edges)
    A, b = assembly(X, Y, etov, lam1=lam1, lam2=lam2, qt=qt, M=M)
    b_neumann1 = neumann_boundary_conditions(X, Y, lam1, lam2, etov, upper_edges, qt, b, test_case=test_case)
    A_final, b_final = dirichlet_boundary_conditions(A, b_neumann1, lower_edges, etov, X, Y, "2_7",
                                                     test_case)
    u_hat = np.linalg.solve(A_final, b_final)
    u = u_true(X, Y, test_case)
    print(list(zip(u_hat, u)))
    print(u_hat)


if __name__ == '__main__':
    main()
