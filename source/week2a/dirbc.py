# function [A,b] = dirbc(bnodes,f,A,b)
#     M = length(b);
#     nodes = setdiff(linspace(1,M,M),bnodes);
#     Nodes = zeros(M,1);
#     Nodes(nodes)=1;
#     for i = bnodes
#         A(i,i) = 1;
#         b(i) = f(i);
#         for j = 1:M
#             if i~=j
#                 A(i,j) = 0;
#                 if Nodes(j)
#                     b(j) = b(j)-A(j,i)*f(i);
#                     A(j,i) = 0;
from collections import defaultdict

import numpy as np


def dirbc(bnodes, f, A1, b1):
    # bnodes or bedges? check!
    A = np.array(A1)
    b = np.array(b1)
    M = len(b)
    for i in bnodes:
        A[i - 1, i - 1] = 1
        b[i - 1] = f[i - 1]

        for j in range(1, M):
            if i != j:
                A[i - 1, j - 1] = 0
                if j - 1 not in bnodes:
                    b[j - 1] = b[j - 1] - A[j - 1, i - 1] * f[i - 1]
                    A[j - 1, i - 1] = 0
    return A, b
