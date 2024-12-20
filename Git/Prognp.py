import numpy as np

def householder_qr(A):
    """QR-разложение методом Хаусхолдера."""
    m, n = A.shape
    Q = np.eye(m)
    R = A.copy()

    for k in range(n):
        # Выделяем столбец
        x = R[k:, k]

        # Формируем вектор отражения
        e1 = np.zeros_like(x)
        e1[0] = np.linalg.norm(x) * (-1 if x[0] < 0 else 1)
        v = x + e1
        v /= np.linalg.norm(v)

        # Формируем матрицу Хаусхолдера
        H_k = np.eye(m)
        H_k[k:, k:] -= 2.0 * np.outer(v, v)

        # Применяем H_k к R и накапливаем Q
        R = H_k @ R
        Q = Q @ H_k.T

    return Q, R

def qr_algorithm(A, max_iter=1000, tol=1e-10):
    """
    QR-алгоритм для нахождения собственных значений.
    
    Аргументы:
    A -- начальная матрица (должна быть квадратной)
    max_iter -- максимальное количество итераций
    tol -- допустимая погрешность для сходимости

    Возвращает:
    A -- верхнетреугольная матрица (приближение формы Шура)
    """
    n, m = A.shape
    assert n == m, "Матрица должна быть квадратной!"

    Ak = A.copy()

    for i in range(max_iter):
        # Шаг 1: QR-разложение методом Хаусхолдера
        Q, R = householder_qr(Ak)

        # Шаг 2: Обновление матрицы
        Ak = R @ Q

        # Шаг 3: Проверка сходимости (мелкие поддиагональные элементы)
        off_diag = np.sqrt(np.sum(np.tril(Ak, -1)**2))
        if off_diag < tol:
            break

    return Ak

# Пример использования
A = np.array([[4.91, -5.0, -14.91],[5.01, -5.1, -14.81],[5.21, -5.3, -14.610000000000001]], dtype=float)

# Применяем QR-алгоритм
result = qr_algorithm(A)

print("Результирующая матрица:")
print(result)
