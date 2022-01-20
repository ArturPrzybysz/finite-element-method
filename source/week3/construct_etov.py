from collections import defaultdict


def construct_element_table(element_count1: int, element_count2: int):
    etov_dict = dict()
    element_to_base = dict()
    base_to_elements = defaultdict(list)
    for e1 in range(element_count1):
        for e2 in range(element_count2):
            e = 1 + 2 * e1 + e2 * (2 * element_count1)
            upper_left = 1 + e1 + e2 * (element_count1 + 1)  #: lewy górny?
            lower_left = upper_left + 1
            upper_right = upper_left + element_count1 + 1
            lower_right = upper_right + 1

            etov_dict[e] = (upper_right, upper_left, lower_right)
            etov_dict[e + 1] = (lower_left, lower_right, upper_left)

            edge = (lower_right, upper_left)
            element_to_base[e] = edge
            element_to_base[e + 1] = edge
            base_to_elements[edge].append(e)
            base_to_elements[edge].append(e + 1)

    M = 2 * (element_count1 * element_count2)
    sorted_etov = {e: etov_dict[e] for e in sorted(etov_dict.keys())}
    return sorted_etov, element_to_base, base_to_elements, M
