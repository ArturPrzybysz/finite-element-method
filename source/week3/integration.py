def integration(A, VX, VY, VN):
    '''
    A_parameter: list with all A
    VX:
    VY
    VN
    '''

    xa = VX[0]
    ya = VY[0]
    xn = VX[1]
    yn = VY[1]
    x = [] #all our domain
    u = u(xa, ya)

    A1 = A[0]
    A2 = A[1]
    A3 = A[2]
    A4 = A[3]
    A5 = A[4]
    A6 = A[5]
    A7 = A[6]
    A8 = A[7]

    integrand = (A2 - A6) ** 2 * (xn - xa) * ((ya + x * (yn - ya) / (xn - xa)) ** 3 - ya ** 3) / 3 + (
                (A2 - A6) * (A1 - A5) * (-xa ** 2 + xn ** 2) + 2 * (A3 * u - A7 * u + A4 - A8) * (A2 - A6) * (
                    xn - xa)) * ((ya + x * (yn - ya) / (xn - xa)) ** 2 - ya ** 2) / 2 + (A1 - A5) ** 2 * (
                     -xa ** 3 + xn ** 3) * x * (yn - ya) / (xn - xa) / 3 + (A3 * u - A7 * u + A4 - A8) * (A1 - A5) * (
                     -xa ** 2 + xn ** 2) * x * (yn - ya) / (xn - xa) + (A3 * u - A7 * u + A4 - A8) ** 2 * x * (yn - ya)

    return integrand