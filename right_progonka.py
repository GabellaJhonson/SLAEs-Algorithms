import numpy as np
eps = 0.00000001
def leftSweepMet(A, b):
    n = len(A)
    m = (n + 1) // 2
    diagonalC = np.zeros(n)
    diagonalA = np.zeros(n)
    diagonalB = np.zeros(n)

    for i in range(n):
        diagonalC[i] = A[i][i]

    for i in range(n - 1):
        diagonalA[i + 1] = -A[i + 1][i]
        diagonalB[i] = -A[i][i + 1]

    print("\nДиагональ A:\n\t", diagonalA)
    print("\nДиагональ B:\n\t", diagonalB)
    print("\nДиагональ C:\n\t", diagonalC)
    alpha = np.zeros(n)
    beta = np.zeros(n)

    alpha[0] = diagonalB[0] / diagonalC[0]
    beta[0] = b[0] / diagonalC[0]

    for i in range(1, n):
        alpha[i] = diagonalB[i] / (diagonalC[i] - diagonalA[i] * alpha[i-1])
        beta[i] = (b[i] + diagonalA[i] * beta[i - 1]) / (diagonalC[i] - diagonalA[i] * alpha[i-1])
        alpha[i] = alpha[i] if abs(alpha[i]) > eps else 0.0
        beta[i] = beta[i] if abs(beta[i]) > eps else 0.0
    print("\nКоэффициенты Alpha:\n\t", alpha)
    print("\nКоэффициенты Beta:\n\t", beta)

    epsilont = np.zeros(n)
    nue = np.zeros(n)

    epsilont[n - 1] = diagonalA[n - 1] / diagonalC[n - 1]
    nue[n - 1] = b[n - 1] / diagonalC[n - 1]

    for i in range(n - 2, m - 1, -1):
        epsilont[i] = diagonalA[i] / (diagonalC[i] - diagonalB[i] * epsilont[i + 1])
        nue[i] = (nue[i + 1] * diagonalB[i] + b[i]) / (diagonalC[i] - diagonalB[i] * epsilont[i + 1])
        epsilont[i] = epsilont[i] if abs(epsilont[i]) > eps else 0.0
        nue[i] = nue[i] if abs(nue[i]) > eps else 0.0

    print("\nКоэффициенты Epsilont:\n\t", epsilont)
    print("\nКоэффициенты Nue:\n\t", nue)

    massX = np.zeros(n)
    massX[m - 1] = ((b[m - 1] + diagonalB[m - 1]*nue[m] + diagonalA[m - 1]*beta[m - 2])/
    (diagonalC[m - 1] - diagonalA[m - 1]*alpha[m - 2] - diagonalB[m - 1]*epsilont[m]))
    
    for i in range(1, m):
        massX[m - i - 1] = alpha[m - i - 1] * massX[m - i] + beta[m - i - 1]
        massX[m + i - 1] = epsilont[m + i - 1] * massX[m + i - 2] + nue[m + i - 1]
    print("\nВектор решений X:\n\t", massX)

    return massX
matrix_A = np.array([[0.6444, 0, 0, 0, 0],
              [-0.0395, 0.4208, 0, 0, 0],
              [0, -0.1184, 0.7627, 0.0145, 0],
              [0, 0, -0.0960, 0.7627, 0],
              [0, 0, 0, -0.0158, 0.5523]])
f = np.array([1.2677, 1.6819, -2.3657, -6.5369, 2.8351])
print(leftSweepMet(matrix_A,f))