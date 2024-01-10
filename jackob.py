import numpy as np
from math import log
epsilon = 10**(-5)

def Jacobi_method(A, b):
    B = np.zeros((5, 5))
    counter = 0
    for i in range(5):
        for j in range(5):
            if i != j:
                B[i][j] = -A[i][j]/A[i][i]
    print(B)
    g = np.zeros(5).transpose()
    for i in range(5):
        g[i] = b[i]/A[i][i]
    print(g)

    x_k = b
    x_k_1 = b.dot(1.2)

    while(np.linalg.norm(x_k_1 - x_k, ord=1) > epsilon) :
        print(f"\n{np.linalg.norm(x_k_1 - x_k, ord=1)} > {epsilon}\n")
        print(f"------------------iteration {counter}------------------\n")

        x_k  = x_k_1
        print(x_k)
        x_k_1 = B.dot(x_k) +g
        print(x_k_1)
        #if counter == 2 : break
        counter+=1

    print(f"\n-------------- Result ------------------\n")
    result = A.dot(x_k_1) - b
    k = log(epsilon*(1-np.linalg.norm(B))/np.linalg.norm(g, 1), np.linalg.norm(B))-1
    print(f"Решение х:{x_k_1}")
    print(f"Количество итераций:{counter}")
    print(f"Расчетное количество итераций:{k}")
    print(f"Вектор невязки:{result}")
    print(f"Норма невязки:{np.linalg.norm(result)}")
    print(f"Норма матрицы B:{np.linalg.norm(B)}")
    
    return B
A = np.array([[0.6444, 0.0000, -0.1683, 0.1184, 0.1973],
              [-0.0395, 0.4208, 0.0000, -0.0802, 0.0263],
              [0.0132, -0.1184, 0.7627, 0.0145, 0.0460],
              [0.0395, 0.0000, -0.0960, 0.7627, 0.0000],
              [0.0263, -0.0395, 0.1907, -0.0158, 0.5523]])
b = np.array([1.2677, 1.6819, -2.3657, -6.5369, 2.8351])
Jacobi_method(A, b)