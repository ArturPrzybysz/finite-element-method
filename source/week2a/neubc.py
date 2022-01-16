from collections import defaultdict

import numpy as np

# function b = neubc(VX,VY,EToV,beds,q,b)
#     for p = 1:length(beds(:,1))
#
#         %n is the number of element to which the edge belongs
#         %r is the local node number of first node on the edge
#         n = beds(p,1);
#         r = beds(p,2);
#
#         % s is the local node number of the second node on the edge
#         % r=1 ⇒ s=2, r=2 ⇒ s=3, r=3 ⇒ s=1
#         if r == 1 || r == 2
#             s = r+1;
#         else
#             s = 1;
#         end
#
#         %i and j are the global node numbers
#         i = EToV(n,r);
#         j = EToV(n,s);
#
#         % Finding the normal vector to calculate the Neuman boundary
#         % λ1⋅ux⋅n1 + λ2⋅uy⋅n2 = −q(x,y)
#         [n1, n2] = outernormal(n,r,VX,VY,EToV);
#
#         %Compute q1(p) and q2(p) from (2.41)
#         if isnumeric(q) == 1
#             q12 = (q/2)*sqrt((VX(j) - VX(i))^2+(VY(j) - VY(i))^2);
#         else
#             xc = 0.5 * (VX(i) + VX(j));
#             yc = 0.5 * (VY(i) + VY(j));
#             q12 = (q(xc,yc,n1,n2)/2)*sqrt((VX(j) - VX(i))^2+(VY(j) - VY(i))^2);
#         end
#
#         % Applying boundary functions onto b
#         b(i) = b(i) - q12;
#         b(j) = b(j) - q12;
#     end
# end
from source.week2a.outernormal import outernormal


def neumann_boundary_conditions(VX, VY, lam1, lam2, EToV, boundary_edges, q, b1, b_target):
    history = defaultdict(list)
    for idx, b in enumerate(b1):
        history[idx].append(b)

    b = np.array(b1)
    for e, r in boundary_edges:
        s = (r + 1) % 3
        vertice_tuple = EToV[e]
        # i and j are the global node numbers
        i = vertice_tuple[r]
        j = vertice_tuple[s]
        #  Finding the normal vector to calculate the Neuman boundary
        #  λ1⋅ux⋅n1 + λ2⋅uy⋅n2 = −q(x,y)
        #    [n1, n2] = outernormal(n,r,VX,VY,EToV);
        [n1, n2] = outernormal(e, r + 1, VX, VY, EToV)

        # Compute q1(p) and q2(p) from (2.41)
        if not hasattr(q, '__call__') == 1:
            q12 = ...  # (q / 2) * np.sqrt((VX(j) - VX(i)) ** 2 + (VY(j) - VY(i)) ** 2)
        else:
            x_j = VX[j - 1]
            x_i = VX[i - 1]
            y_j = VY[j - 1]
            y_i = VY[i - 1]
            xc = 0.5 * (VX[i - 1] + VX[j - 1])
            yc = 0.5 * (VY[i - 1] + VY[j - 1])
            q12 = (q(xc, yc, n1, n2) / 2) * np.sqrt((x_j - x_i) ** 2 + (y_j - y_i) ** 2)

        b[i - 1] = b[i - 1] - q12
        history[i-1].append(f"- {q12}")
        target1 = b_target[i - 1]
        b[j - 1] = b[j - 1] - q12
        history[j-1].append(f"- {q12}")

        target2 = b_target[j - 1]

    return b
