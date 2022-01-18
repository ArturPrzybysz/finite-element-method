def integration(A, VX, VY, VN):
    '''
    A_parameter: list with all A
    VX:
    VY
    VN
    '''

    XA = VX[0]
    YA = VY[0]
    XN = VX[1]
    YN = VY[1]
    X = [] #all our domain

    B1 = A[0] - A[3] #A1 - A4
    B2 = A[1] - A[4]  # A2 - A5
    B3 = A[2] - A[5]  # A3 - A6

    integrand_part1 = (((B1*XA)+(B2*XA*(YN*XA-YA*XA))/(XN*XA-XA**2))*(X*XA)+B2*XA*YA*XA+B3*XA)**4

    integrand_part2 = 12*B2*XA*(B1*XN+())

    return