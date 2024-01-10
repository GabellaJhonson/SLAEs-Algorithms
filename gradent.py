import numpy as np

def gradient_descent(A, b, tolerance):
    n = len(b)
    E = np.identity(n)
    At = A.transpose() # находим Аt
    A = np.dot(At, A) # перемножаем А на Аt, теперь А - симметрическая
    b = np.dot(At, b) # то же самое с b
    xk = b # начальное приближение
    x = np.zeros(n)
    k = 0
    while True: # итерационный процесс
        rk = np.dot(A, xk) - b # считаем rk
        x = xk - np.dot(rk, np.dot(rk, rk) / np.dot(np.dot(A, rk), rk))
        k += 1
        if abs(np.linalg.norm(x, 1) - np.linalg.norm(xk, 1)) < tolerance: # условие окончания итерации
            break
        xk = x
    r = np.dot(A, x) - b # находим вектор невязки
    rnorm = np.linalg.norm(r, 1)
    print(f"Решение х:{x}")
    print(f"Количество итераций:{k}")
    print(f"Вектор невязки:{r}")
    print(f"Норма невязки:{rnorm}")

# Пример использования
A = np.array([[0.6444, 0.0000, -0.1683, 0.1184, 0.1973],
              [-0.0395, 0.4208, 0.0000, -0.0802, 0.0263],
              [0.0132, -0.1184, 0.7627, 0.0145, 0.0460],
              [0.0395, 0.0000, -0.0960, 0.7627, 0.0000],
              [0.0263, -0.0395, 0.1907, -0.0158, 0.5523]])
b = np.array([1.2677, 1.6819, -2.3657, -6.5369, 2.8351])
tolerance = 1e-5  # Предел точности
# Шаг обучения = 1

gradient_descent(A, b, tolerance)
#[ 0.99821505  1.99986528 -2.99975971 -9.00000843  6.00705353]