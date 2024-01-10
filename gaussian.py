import math
import numpy as np
from decimal import Decimal
def gaussian(matrix_, inhomogeneity):
    matrix = matrix_.copy() # Копируем исходную матрицу в локальную переменную
    insertions = 0 # Число перестановок
    indexes = [i for i in range (matrix.shape[0])] # Вектор индексов
    unit_matrix = np.eye(matrix.shape[0]) # Единичная матрица
    
    # Выбор главного элемента по строке
    for k in range(matrix.shape[0]):
        leading_column = k # Столбец, в котором ведущий элемент, считаем равным номеру шага
        for i in range(k, matrix.shape[0]): # находим максимальный по модулю элемент в строке
            if math.fabs(matrix[k][i]) > math.fabs(matrix[leading_column][k]): # Если нашли в строке элемент ...
                # ...больший, чем ведущий, то..
                leading_column = i # Номер столбца с ведущим элементом принимает значение того, ...
                # ...в котором находится больший элемент, чем текущий ведущий элемент
        for j in range(k, matrix.shape[1]): # меняем местами столбец, в котором главный элемент, со столбцом равном номеру шага
            matrix[j][leading_column], matrix[j][k] =  matrix[j][k], matrix[j][leading_column]
        insertions+=1 # повышаем число перестановок
        indexes[leading_column], indexes[k] = indexes[k], indexes[leading_column]
        
        # Прямой ход (полностью дублирует формулу прямого хода, указанную выше)
        q = inhomogeneity[k] / matrix[k][k] 
        for j in range(matrix.shape[0] - 1, k - 1, -1):
            c = matrix[k][j] / matrix[k][k]
            for i in range(matrix.shape[0] - 1, k, -1):
                matrix[i][j] = matrix[i][j] - matrix[i][k]*c
                if j == matrix.shape[0] - 1:
                    inhomogeneity[i] = inhomogeneity[i] - matrix[i][k]*q
                
    
    # Обратный ход (полностью дублирует формулу обратного хода, указанную выше)
    results_with_insertions = np.zeros(matrix.shape[0]) # Нулевой вектор
    for i in range(matrix.shape[0]-1, -1, -1):
        summary = 0
        for j in range(i+1, matrix.shape[0]):
            summary += matrix[i][j]*results_with_insertions[j]
        results_with_insertions[i] = (inhomogeneity[i] - summary) / matrix[i][i]
        
    
    # Перенумерация индексов в векторе решений
    results = np.zeros(matrix.shape[0]) # Создаем итоговый вектор решений, заполненный нулями
    for i in range(matrix.shape[0]):
        results[indexes[i]] = results_with_insertions[i] # Переносим значения из получившегося столбца решений ...
        # в созданный в той последовательности, в которой элементы столбца должны расоплагаться
        
    #Вычисление вектора невязки
    discrepancy_vector = np.zeros(matrix.shape[0]) # Создаем вектор заполненный нулями
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            discrepancy_vector[i] += matrix[i][j]*results[j]; # r += Ax
        discrepancy_vector[i] -= inhomogeneity[i] # r -= b
    
    # Вычисление определителя
    determinant = (-1)**insertions # Создаем переменную, в которой содержится (-1) в степени числа перестановок
    for k in range(matrix.shape[0]):
        determinant *= matrix[k][k] # Перемножаем все диагональные элементы получившейся треугольной матрицы
        
    return results, discrepancy_vector, determinant; # Возвращаем значения решения, невязки и определителя
def inverse_A(A):
    n = 5
    E = np.eye(n)
    A_inv = np.empty((0, 5))
    for i in range(n):
        A_inv = np.vstack((A_inv, gaussian(A, E[i])[0]))
        # A_inv.append(gaussian(A, E[i])[0])
    return A_inv.T
    # print(*A_inv, sep='\n')
def condition_N(A):
    nu = np.linalg.norm(A, 2) * np.linalg.norm(inverse_A(A), 2)
    return nu  
A = np.array([[0.6444, 0.0000, -0.1683, 0.1184, 0.1973],
              [-0.0395, 0.4208, 0.0000, -0.0802, 0.0263],
              [0.0132, -0.1184, 0.7627, 0.0145, 0.0460],
              [0.0395, 0.0000, -0.0960, 0.7627, 0.0000],
              [0.0263, -0.0395, 0.1907, -0.0158, 0.5523]])
b = np.array([1.2677, 1.6819, -2.3657, -6.5369, 2.8351])
# решение
print("solution \n",gaussian(A,b)[0])
# решение
print("discrepancy vector \n",gaussian(A,b)[1])
# решение
print("determinat \n ",gaussian(A,b)[2])
# обратная
print("A^-1 \n",inverse_A(A))
# матрица невзяки A*A^-1 - E
print("A*A^-1 - E \n",np.dot(inverse_A(A), A) - np.eye(5))
# число обусловленности
print("conditional number ",condition_N(A))
matrix_A = np.array([[0.6444, 0, 0, 0, 0],
              [-0.0395, 0.4208, 0, 0, 0],
              [0, -0.1184, 0.7627, 0.0145, 0],
              [0, 0, -0.0960, 0.7627, 0],
              [0, 0, 0, -0.0158, 0.5523]])
f = np.array([1.2677, 1.6819, -2.3657, -6.5369, 2.8351])
# решение
print("solution \n",gaussian(matrix_A,f)[0])
# решение
print("discrepancy vector \n",gaussian(matrix_A,f)[1])
decimal_x = [Decimal(str(i)) for i in gaussian(matrix_A,f)[0]]
decimal_f = [Decimal(str(i)) for i in f]
decimal_matrix = np.vectorize(lambda i: Decimal(str(i)))(matrix_A)
# np.set_printoptions(precision=22, suppress=True)
print("\nВектор невязки: \n\t", gaussian(matrix_A, f)[1])
print("\Норма невязки: \n\t",np.linalg.norm(gaussian(matrix_A, f)[1]))
# print("\nВектор невязки: \n\t", gaussian(decimal_matrix, decimal_f)[1])
# print("\Норма невязки: \n\t",np.linalg.norm(gaussian(decimal_matrix, decimal_f)[1]))