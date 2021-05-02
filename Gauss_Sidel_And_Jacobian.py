import sys

def determinant_recursive(A, total=0):
    # Section 1: store indices in list for row referencing
    indices = list(range(len(A)))

    # Section 2: when at 2x2 submatrices recursive calls end
    if len(A) == 2 and len(A[0]) == 2:
        val = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return val

    # Section 3: define submatrix for focus column and
    #      call this function
    for fc in indices:  # A) for each focus column, ...
        # find the submatrix ...
        As = A[:]  # B) make a copy, and ...
        As = As[1:]  # ... C) remove the first row
        height = len(As)  # D)

        for i in range(height):
            # E) for each remaining row of submatrix ...
            #     remove the focus column elements
            As[i] = As[i][0:fc] + As[i][fc + 1:]

        '''
             -  +
        +[3, 4, 7]
             +  -        
        -[8, 4, 8]
             -  +
        +[2, 1, 1]
        '''

        sign = (-1) ** (fc % 2)  # F)
        # G) pass submatrix recursively
        sub_det = determinant_recursive(As)
        # H) total all returns from recursion
        total += sign * A[0][fc] * sub_det

    return total

def zeros_matrix(rows, cols):
    A = []
    for i in range(rows):
        A.append([])
        for j in range(cols):
            A[-1].append(0.0)

    return A

def matrix_multiply(A, B):
    rowsA = len(A)
    colsA = len(A[0])

    rowsB = len(B)
    colsB = len(B[0])

    if colsA != rowsB:
        print('Number of A columns must equal number of B rows.')
        sys.exit()

    C = zeros_matrix(rowsA, colsB)

    for i in range(rowsA):
        for j in range(colsB):
            total = 0
            for ii in range(colsA):
                total += A[i][ii] * B[ii][j]
            C[i][j] = total

    return C

def copy_matrix(M):
    rows = len(M)
    cols = len(M[0])

    MC = zeros_matrix(rows, cols)

    for i in range(rows):
        for j in range(rows):
            MC[i][j] = M[i][j]

    return MC

def identity_matrix(n):
    I = []
    for i in range(n):
        I.append([])
        for j in range(n):
            if i == j:
                I[i].append(1)
            else:
                I[i].append(0)
    return I

def invert_matrix(A, tol=None):
    """
    Returns the inverse of the passed in matrix.
        :param A: The matrix to be inversed

        :return: The inverse of the matrix A
    """
    # Section 1: Make sure A can be inverted.

    # Section 2: Make copies of A & I, AM & IM, to use for row ops
    n = len(A)
    AM = copy_matrix(A)
    I = identity_matrix(n)
    IM = copy_matrix(I)

    # Section 3: Perform row operations
    indices = list(range(n))  # to allow flexible row referencing ***
    for fd in range(n):  # fd stands for focus diagonal
        fdScaler = 1.0 / AM[fd][fd]
        # FIRST: scale fd row with fd inverse.
        for j in range(n):  # Use j to indicate column looping.
            AM[fd][j] *= fdScaler
            IM[fd][j] *= fdScaler
        # SECOND: operate on all rows except fd row as follows:
        for i in indices[0:fd] + indices[fd + 1:]:
            # *** skip row with fd in it.
            crScaler = AM[i][fd]  # cr stands for "current row".
            for j in range(n):
                # cr - crScaler * fdRow, but one element at a time.
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
                IM[i][j] = IM[i][j] - crScaler * IM[fd][j]

    # Section 4: Make sure IM is an inverse of A with specified tolerance
    # if check_matrix_equality(I, matrix_multiply(A, IM), tol):
    return IM
    # else:
    #     raise ArithmeticError("Matrix inverse out of tolerance.")

def print_matrix(Title, M):
    print(Title)
    for row in M:
        print([round(x, 3) + 0 for x in row])

def norma(A):
    B = A[:]
    B = [[abs(element) for element in row] for row in B]
    B = list(map(lambda row: sum(row), B))
    return max(B)

def isDominant(A):
    for i in range(len(A)):
        dominant = A[i][i]
        s = sum(A[i]) - dominant
        if s >= dominant:
            return False
    return True

def getLDU(A):
    L = []
    D = []
    U = []
    for i in range(len(A)):
        L.append([])
        D.append([])
        U.append([])
        for j in range(len(A[i])):
            if i > j:
                L[i].append(A[i][j])
                D[i].append(0)
                U[i].append(0)
            elif i == j:
                L[i].append(0)
                D[i].append(A[i][j])
                U[i].append(0)
            else:
                L[i].append(0)
                D[i].append(0)
                U[i].append(A[i][j])
    return L, D, U

def Jaccoian(A, b, epsilon):
    dominant = isDominant(A)
    size = len(A)
    if(dominant):
        print(f'Since the given matrix has dominant diagonal Jaccoian will converge.')
    else:
        print(f'Since the given matrix does not has dominant diagonal Jaccoian might diverge.')

    L, D, U = getLDU(A)
    D_inverse = invert_matrix(D)
    print_matrix("A", A)
    print_matrix("L", L)
    print_matrix("D", D)
    print_matrix("U", U)
    print_matrix("D inverse", D_inverse)
    minus_D_inverse = [[-D_inverse[i][j] for i in range(size)] for j in range(size)]
    print_matrix("-D^-1", minus_D_inverse)
    L_plus_U = []
    for i in range(size):
        L_plus_U.append([])
        for j in range(len(L)):
            if i != j:
                L_plus_U[i].append(L[i][j] + U[i][j])
            else:
                L_plus_U[i].append(0)

    print_matrix("L + U", L_plus_U)
    G = matrix_multiply(minus_D_inverse, L_plus_U)
    print_matrix("G", G)

    cond = norma(G)
    if cond >= 1:
        raise Exception(f'The norm of the matrix is {cond} which is greater or equal to 1\nTherefore gauss sidel diverge.')
    xr = [0] * len(A)

    def f(A, i, *args):
        ans = b[i][0]
        for j in range(len(A)):
            if i != j:
                ans -= A[i][j] * args[0][j]
        ans /= A[i][i]
        return ans

    xr_plus_1 = [2 * epsilon] * size
    tmp = xr[:]
    iteration = 1
    while abs(xr_plus_1[0] - tmp[0]) > epsilon:
        # print xr
        tmp = xr[:]
        print(f'Iteration {iteration}')
        print("xr:")
        for i in range(size):
            print(f'x{i}', end='\t\t\t\t\t')
        print("")
        for i in range(size):
            print("{:<20} ".format(xr[i]),end='')

        # calculate xr+1
        for i in range(size):
            xr_plus_1[i] = f(A, i, xr)
        xr = xr_plus_1[:]
        # print("xr+1:")
        print("")
        print("xr+1:")
        for i in range(size):
            print(f'x{i}', end='\t\t\t\t\t')
        print("")
        for i in range(size):
            print("{:<20} ".format(xr_plus_1[i]),end='')
        print("")
        iteration += 1
    if not (dominant):
        print("Although the matrix does not have dominant diagonal it converge to\n{0}".format(xr_plus_1))

def gaussSidel(A, b, epsilon):
    dominant = isDominant(A)
    size = len(A)
    if(dominant):
        print(f'Since the given matrix has dominant diagonal gauss sidel will converge.')
    else:
        print(f'Since the given matrix does not has dominant diagonal gauss sidel might diverge.')

    L, D, U = getLDU(A)
    print_matrix("A", A)
    print_matrix("L", L)
    print_matrix("D", D)
    print_matrix("U", U)
    L_plus_D = []

    for i in range(size):
        L_plus_D.append([])
        for j in range(len(L)):
            if i == j:
                L_plus_D[i].append(D[i][j])
            elif i > j:
                L_plus_D[i].append(L[i][j])
            else:
                L_plus_D[i].append(0)

    print_matrix("L + D", L_plus_D)
    L_plus_D_inverse = invert_matrix(L_plus_D)
    minus_L_plus_D_inverse = [[-L_plus_D_inverse[i][j] for i in range(size)] for j in range(size)]

    G = matrix_multiply(minus_L_plus_D_inverse, U)
    print_matrix("G", G)

    cond = norma(G)
    if cond >= 1:
        raise Exception(f'The norm of the matrix is {cond} which is greater or equal to 1\nTherefore gauss sidel diverge.')
    xr = [0] * len(A)

    def f(A, i, *args):
        ans = b[i][0]
        for j in range(len(A)):
            if i != j:
                ans -= A[i][j] * args[0][j]
        ans /= A[i][i]
        return ans

    xr_plus_1 = [2 * epsilon] * size

    tmp = xr[:]
    iteration = 1

    while abs(xr_plus_1[0] - tmp[0]) > epsilon:
        tmp = xr[:]
        # print xr
        print(f'Iteration {iteration}')
        print("xr:")
        for i in range(size):
            print(f'x{i}', end='\t\t\t\t\t')
        print("")
        for i in range(size):
            print("{:<20} ".format(xr[i]),end='')

        # calculate xr+1
        for i in range(size):
            xr_plus_1[i] = f(A, i, xr)
            xr[i] = xr_plus_1[i]
        # print("xr+1:")
        print("")
        print("xr+1:")
        for i in range(size):
            print(f'x{i}', end='\t\t\t\t\t')
        print("")
        for i in range(size):
            print("{:<20} ".format(xr_plus_1[i]),end='')
        print("")
        iteration += 1
    if not(dominant):
        print("Although the matrix does not have dominant diagonal it converge to\n{0}".format(xr_plus_1))



A = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]
b = [[2], [6], [5]]
epsilon = 0.001

choice = -1
while choice not in [1, 2]:
    print(f'A = {A}')
    print(f'b = {b}')
    choice = int(input("Please enter one of the choices:\n1) Gauss Sidel\n2) Jacobian"))
    if choice == 1:
        gaussSidel(A, b, epsilon)
    elif choice == 2:
        Jaccoian(A, b, epsilon)
    else:
        continue



