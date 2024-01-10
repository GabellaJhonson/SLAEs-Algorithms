import numpy as np

def leftSweepMet(A, b):
    n = len(A)
    
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

    alpha[n - 1] = diagonalA[n - 1] / diagonalC[n - 1]
    beta[n - 1] = b[n - 1] / diagonalC[n - 1]

    for i in range(n - 2, 0, -1):
        alpha[i] = diagonalA[i] / (diagonalC[i] - diagonalB[i] * alpha[i + 1])
        beta[i] = (beta[i + 1] * diagonalB[i] + b[i]) / (diagonalC[i] - diagonalB[i] * alpha[i + 1])

    print("\nКоэффициенты Alpha:\n\t", alpha)
    print("\nКоэффициенты Beta:\n\t", beta)

    massX = np.zeros(n)
    massX[0] = (b[0] + beta[1] * diagonalB[0]) / (diagonalC[0] - alpha[1] * diagonalB[0])

    for i in range(1, n):
        massX[i] = alpha[i] * massX[i - 1] + beta[i]

    print("\nВектор решений X:\n\t", massX)

    return massX
matrix_A = np.array([[0.6444, 0, 0, 0, 0],
              [-0.0395, 0.4208, 0, 0, 0],
              [0, -0.1184, 0.7627, 0.0145, 0],
              [0, 0, -0.0960, 0.7627, 0],
              [0, 0, 0, -0.0158, 0.5523]])
f = np.array([1.2677, 1.6819, -2.3657, -6.5369, 2.8351])
print(leftSweepMet(matrix_A,f))