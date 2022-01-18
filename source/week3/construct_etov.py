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
    M = (element_count1 + 1) * (element_count2 + 1)
    return etov_dict, M
