from typing import Tuple


def part_of_integral(plane1, plane2, A, N, x):
    A1, A2, A3, A4 = plane1
    A5, A6, A7, A8 = plane2
    xn, yn = N
    xa, ya = A
    return (-A2 / A3 + A6 / A7) ** 2 * (xn - xa) * ((ya + x * (yn - ya) / (xn - xa)) ** 3 - ya ** 3) / 3 + (
            (-A2 / A3 + A6 / A7) * (-A1 / A3 + A5 / A7) * (-xa ** 2 + xn ** 2) + 2 * (-A4 / A3 + A8 / A7) * (
            -A2 / A3 + A6 / A7) * (xn - xa)) * ((ya + x * (yn - ya) / (xn - xa)) ** 2 - ya ** 2) / 2 + (
                   -A1 / A3 + A5 / A7) ** 2 * (-xa ** 3 + xn ** 3) * x * (yn - ya) / (xn - xa) / 3 + (
                   -A4 / A3 + A8 / A7) * (-A1 / A3 + A5 / A7) * (-xa ** 2 + xn ** 2) * x * (yn - ya) / (xn - xa) + (
                   -A4 / A3 + A8 / A7) ** 2 * x * (yn - ya)


def integration(plane1, plane2, A, N, integration_range: Tuple):
    int1 = part_of_integral(plane1, plane2, A, N, integration_range[1])
    int2 = part_of_integral(plane1, plane2, A, N, integration_range[0])
    result = int1 - int2
    return result
