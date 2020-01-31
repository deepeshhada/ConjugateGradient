from copy import deepcopy


def GetAlpha(rk, pk, A):
    return innerProduct(rk, rk) / innerProduct(pk, MatrixVectorProduct(A, pk))


def GetXkp1(xk, pk, alpha):
    return AddVectors(xk, ScalarVectorProduct(alpha, pk))


def GetRkp1(rk, pk, A, alpha):
    return AddVectors(rk, ScalarVectorProduct(alpha, MatrixVectorProduct(A, pk)))


def GetBeta(rk, rkp1):
    print("GetBeta ", innerProduct(rkp1, rkp1))
    return innerProduct(rkp1, rkp1) / innerProduct(rk, rk)


def GetPkp1(pk, rkp1, beta):
    return SubtractVectors(ScalarVectorProduct(beta, pk), rkp1)


def VectorNorm(vector):
    result = 0
    for i in range(len(vector)):
        result += vector[i] ** 2
    return (result) ** 0.5


def AddVectors(vector1, vector2):
    vector = [vector1[i] + vector2[i] for i in range(len(vector1))]
    return vector


def SubtractVectors(vector1, vector2):
    # print(vector1)
    # print(vector2)
    vector = [vector1[i] - vector2[i] for i in range(len(vector1))]
    return vector


def innerProduct(vector1, vector2):
    result = 0
    for i in range(len(vector1)):
        result += vector1[i] * vector2[i]
    return result


def ScalarVectorProduct(scalar, vector1):
    vector = [scalar * vector1[i] for i in range(len(vector1))]
    return vector


def MatrixVectorProduct(matrix, vector1):
    vector = []
    for i in range(len(vector1)):
        vector.append(innerProduct(matrix[i], vector1))
    return vector


def createHilbertMatrix(dimension):
    matrix = []
    row = []
    for i in range(0, dimension):
        for j in range(0, dimension):
            row.append(1 / (i + j + 1))
        matrix.append(row)
        row = []

    return matrix


def createOneVector(dimension):
    vector = []
    for i in range(0, dimension):
        vector.append(1)
    return vector


def createZeroVector(dimension):
    vector = []
    for i in range(0, dimension):
        vector.append(0)
    return vector


def StartConjugate():
    dimension = 5

    A = createHilbertMatrix(dimension)
    b = createOneVector(dimension)
    x0 = createZeroVector(dimension)

    r0 = SubtractVectors(MatrixVectorProduct(A, x0), b)  # Calculate r0 = Ax0 - b
    k = 0
    beta = 0
    alpha = 0

    rk = r0
    rkp1 = []
    xk = createZeroVector(dimension)
    xkp1 = []
    pk = ScalarVectorProduct(-1, r0)
    pkp1 = []

    errorList = []

    while VectorNorm(rk) > 10 ** -6:
        errorList.append(VectorNorm(rk))
        alpha = GetAlpha(rk, pk, A)
        xkp1 = GetXkp1(xk, pk, alpha)
        rkp1 = GetRkp1(rk, pk, A, alpha)
        beta = GetBeta(rk, rkp1)
        pkp1 = GetPkp1(pk, rkp1, beta)

        k += 1
        xk = xkp1
        rk = rkp1
        pk = pkp1

        print("k: ", k)
        print("rk Norm: ", VectorNorm(rk))
    errorList.append(VectorNorm(rk))
    print(errorList)
    return


StartConjugate()
