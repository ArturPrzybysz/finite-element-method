import numpy as np

"""
(2.49):  Ω = {(x, y) ∈ R 2 |x 0 ≤ x ≤ x 0 + L_1 , y 0 ≤ y ≤ y 0 + L_2 }
"""


def xy(x0, y0, L1, L2, noelms1, noelms2):
    all_nodes = (1 + noelms1) * (1 + noelms2)
    X = np.zeros(all_nodes)

    for i, x in enumerate(np.linspace(start=x0, stop=x0 + L1, num=noelms1 + 1)):
        start_idx = i * (noelms2 + 1)
        stop_idx = start_idx + (noelms2 + 1)
        X[start_idx:stop_idx] = x
    # for j, y in np.linspace(start=y0, stop=x0 + L1, num=noelms1):
    #     X[i:(i*noelms2)] = x

    first_row = np.linspace(start=y0, stop=y0 + L2, num=noelms2 + 1)[::-1]
    Y = np.tile(first_row, noelms1 + 1)
    return X, Y


def construct_element_table(element_count1: int, element_count2: int):
    etov_dict = dict()

    triangle_to_count = {"upper": 1, "lower": 2}

    for e1 in range(element_count1):
        for e2 in range(element_count2):
            for triangle in ["upper", "lower"]:
                count = triangle_to_count[triangle]
                e = count + 2 * e2 + (2 + element_count1) * e1
                if triangle == "upper":
                    v2 = (e1 * (element_count2 + 1)) + e2 + 1
                    v1 = v2 + element_count2 + 1
                    v3 = v1 + 1
                else:
                    v1 = 2 + (e1 * (element_count2 + 1)) + e2
                    v3 = 1 + (e1 * (element_count2 + 1)) + e2
                    v2 = v3 + element_count2 + 2
                etov_dict[e] = (v1, v2, v3)
    return etov_dict


def main():
    """Test case:
    (x_0 , y_0 ) = (−2.5, −4.8), L_1 = 7.6, L_2 = 5.9, noelms1 = 4, noelms2 = 3
    """
    VX, VY = xy(-2.5, -4.8, 7.6, 5.9, 4, 3)
    etov_dict = construct_element_table(4, 3)
    print()


if __name__ == '__main__':
    main()
