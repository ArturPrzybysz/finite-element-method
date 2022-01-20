from source.week3.visualisation import plot_2d_grid


def sorted_tuple(tpl):
    return tuple(sorted(tpl))


def newest_node_bisection(e_k, EToV, X, Y, element_to_base, edge_to_element):
    base_k = element_to_base[e_k]

    if len(edge_to_element[base_k]) == 1:  # the base of e_k is a boundary edge
        EToV, X, Y, element_to_base, edge_to_element = replace_element_with_two(e_k, EToV, X, Y, element_to_base,
                                                                                edge_to_element)
    else:
        e_j = [element for element in edge_to_element[base_k] if element != e_k][0]
        base_j = element_to_base[e_j]
        if base_j == base_k:
            # Replace each of e k and e j by two sub-triangles
            EToV, X, Y, element_to_base, edge_to_element = replace_element_with_two(e_j, EToV, X, Y, element_to_base,
                                                                                    edge_to_element)
            EToV, X, Y, element_to_base, edge_to_element = replace_element_with_two(e_k, EToV, X, Y, element_to_base,
                                                                                    edge_to_element)
        else:
            EToV, X, Y, element_to_base, edge_to_element = \
                newest_node_bisection(e_j, EToV, X, Y, element_to_base, edge_to_element)
            EToV, X, Y, element_to_base, edge_to_element = \
                newest_node_bisection(e_k, EToV, X, Y, element_to_base, edge_to_element)

    return EToV, X, Y, element_to_base, edge_to_element


def replace_element_with_two(element_idx, EToV, X, Y, element_to_base, edge_to_element):
    r, s, t = EToV[element_idx]
    pt1, pt2 = element_to_base[element_idx]

    r -= 1
    s -= 1
    t -= 1

    pt1 -= 1
    pt2 -= 1
    x_r, x_s, x_t, x_1, x_2 = X[r], X[s], X[t], X[pt1], X[pt2]
    y_r, y_s, y_t, y_1, y_2 = Y[r], Y[s], Y[t], Y[pt1], Y[pt2]
    x_c = (x_1 + x_2) / 2
    y_c = (y_1 + y_2) / 2

    X.append(x_c)
    Y.append(y_c)
    c = len(X)

    last_idx = len(EToV)

    triangle1 = (c, t + 1, r + 1)
    triangle2 = (c, r + 1, s + 1)

    element_to_base[element_idx] = sorted_tuple((t + 1, r + 1))
    element_to_base[last_idx + 1] = sorted_tuple((r + 1, s + 1))

    edge_to_element[sorted_tuple((r + 1, s + 1))].remove(element_idx)
    edge_to_element[sorted_tuple((s + 1, t + 1))].remove(element_idx)
    edge_to_element[sorted_tuple((t + 1, r + 1))].remove(element_idx)

    for edge in [(c, t + 1), (t + 1, r + 1), (r + 1, c)]:
        edge_to_element[sorted_tuple(edge)].add(element_idx)

    for edge in [(c, r + 1), (r + 1, s + 1), (s + 1, c)]:
        edge_to_element[sorted_tuple(edge)].add(last_idx + 1)

    EToV[element_idx] = triangle1
    EToV[last_idx + 1] = triangle2

    # plot_2d_grid(X, Y, EToV, elements_to_plot=list(EToV.keys()))

    return EToV, X, Y, element_to_base, edge_to_element
