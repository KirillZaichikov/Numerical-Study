import numpy as np
import math as m
from scipy.optimize import root
from typing import Callable, Optional, Tuple
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True, precision=50)

class PeriodicPointFinder3D:
    """
    Класс для поиска периодических точек трехмерного отображения.
    """
    
    def __init__(self, 
                 mapping: Callable[[np.ndarray], np.ndarray],
                 period: int,
                 eps: float = 1e-9):
        """
        Инициализация поисковика периодических точек.
        
        Parameters:
        -----------
        mapping : Callable
            Трехмерное отображение F: R^3 -> R^3
        period : int
            Период точек для поиска
        eps : float
            Малое число для численного дифференцирования
        """
        self.mapping = mapping
        self.period = period
        self.eps = eps
    
    def numerical_jacobian(self, x: np.ndarray) -> np.ndarray:
        """
        Численное вычисление якобиана 3x3 в точке x.
        
        Parameters:
        -----------
        x : np.ndarray
            Точка в R^3
            
        Returns:
        --------
        J : np.ndarray
            Матрица Якоби 3x3
        """
        n = len(x)
        J = np.zeros((n, n))
        
        for i in range(n):
            # Вектор возмущения вдоль i-ой координаты
            h = np.zeros(n)
            h[i] = self.eps
            
            # Вычисление производной по i-ой переменной
            f_plus = self.mapping(x + h)
            f_minus = self.mapping(x - h)
            J[:, i] = (f_plus - f_minus) / (2 * self.eps)
        
        return J
    
    def period_mapping(self, x0: np.ndarray) -> np.ndarray:
        """
        Применяет отображение period раз.
        
        Parameters:
        -----------
        x0 : np.ndarray
            Начальная точка
            
        Returns:
        --------
        x_period : np.ndarray
            Точка после period итераций
        """
        x = x0.copy()
        for _ in range(self.period):
            x = self.mapping(x)
        return x
    
    def period_jacobian(self, x0: np.ndarray) -> np.ndarray:
        """
        Вычисляет якобиан period-кратной итерации отображения.
        
        Parameters:
        -----------
        x0 : np.ndarray
            Начальная точка
            
        Returns:
        --------
        J_total : np.ndarray
            Матрица Якоби композиции period отображений
        """
        x = x0.copy()
        J_total = np.eye(3)
        
        for _ in range(self.period):
            # Якобиан текущей итерации
            J_step = self.numerical_jacobian(x)
            # Цепное правило
            J_total = J_step @ J_total
            # Переход к следующей точке
            x = self.mapping(x)
        
        return J_total
    
    def find_periodic_point(self,
                           initial_guess: np.ndarray,
                           method: str = 'hybr',
                           tol: float = 1e-12,
                           maxiter: int = 2000) -> dict:
        """
        Находит периодическую точку заданного периода.
        
        Parameters:
        -----------
        initial_guess : np.ndarray
            Начальное приближение
        method : str
            Метод решения уравнений (из scipy.optimize.root)
        tol : float
            Точность решения
        maxiter : int
            Максимальное число итераций
            
        Returns:
        --------
        result : dict
            Результаты поиска:
            - 'success': успех поиска
            - 'point': найденная точка
            - 'stability': собственные числа
            - 'jacobian': матрица Якоби
            - 'iterations': число итераций
        """
        
        def equation(x: np.ndarray) -> np.ndarray:
            """Уравнение для поиска периодической точки: F^period(x) - x = 0"""
            return self.period_mapping(x) - x
        
        # Решение уравнения F^period(x) = x
#         sol = root(equation, initial_guess, method=method, 
#                    tol=tol, options = {
#     'xtol': 1e-2,       # точность по x
#     'maxfev': 1000,     # максимум вызовов функции
#     'diag': None,       # масштабирование
# })
        sol = root(equation, initial_guess, method=method, tol=tol)
        
        result = {
            'success': sol.success,
            'point': sol.x,
            'iterations': sol.nfev,
            'message': sol.message,
            'residual': np.linalg.norm(sol.fun)
        }
        
        if sol.success:
            # Вычисляем матрицу Якоби и её собственные числа
            J = self.period_jacobian(sol.x)
            # eigenvalues = np.linalg.eigvals(J)
            eigenvalues = np.linalg.eigvals(J)
            
            result['jacobian'] = J
            result['eigenvalues'] = eigenvalues
            result['stability'] = self.classify_stability(eigenvalues)
            result['distance'] = np.linalg.norm(initial_guess - sol.x)
            
            # Проверка, что это действительно точка нужного периода
            if not self.verify_period(sol.x):
                result['success'] = False
                result['message'] = 'Found point has lower period than requested'
        
        return result
    
    def find_multiple_points(self,
                            initial_guesses: np.ndarray,
                            **kwargs) -> list:
        """
        Ищет несколько периодических точек из разных начальных приближений.
        
        Parameters:
        -----------
        initial_guesses : np.ndarray
            Массив начальных приближений (N x 3)
        **kwargs : dict
            Дополнительные аргументы для find_periodic_point
            
        Returns:
        --------
        results : list
            Список результатов
        """
        results = []
        for guess in initial_guesses:
            result = self.find_periodic_point(guess, **kwargs)
            results.append(result)
        return results
    
    def classify_stability(self, eigenvalues: np.ndarray) -> str:
        """
        Классифицирует устойчивость периодической точки.
        
        Parameters:
        -----------
        eigenvalues : np.ndarray
            Собственные числа матрицы Якоби
            
        Returns:
        --------
        stability_type : str
            Тип устойчивости
        """
        # Вычисляем модули собственных чисел
        magnitudes = np.abs(eigenvalues)
        
        if np.all(magnitudes < 1):
            return 'stable (sink)'
        elif np.all(magnitudes > 1):
            return 'unstable (source)'
        elif np.any(magnitudes < 1) and np.any(magnitudes > 1):
            return 'saddle'
        else:
            # Есть собственные числа на единичной окружности
            on_circle = np.sum(np.abs(magnitudes - 1) < 1e-8)
            if on_circle == 1:
                return 'non-hyperbolic (one eigenvalue on unit circle)'
            elif on_circle == 2:
                return 'non-hyperbolic (two eigenvalues on unit circle)'
            else:
                return 'non-hyperbolic'
    
    def verify_period(self, point: np.ndarray, tol: float = 1e-3) -> bool:
        """
        Проверяет, что точка действительно имеет заданный период.
        
        Parameters:
        -----------
        point : np.ndarray
            Проверяемая точка
        tol : float
            Допуск
            
        Returns:
        --------
        is_correct_period : bool
            True если точка имеет правильный период
        """
        x = point.copy()
        
        # Проверяем все меньшие периоды
        for p in range(1, self.period):
            x_temp = point.copy()
            for _ in range(p):
                x_temp = self.mapping(x_temp)
            
            if np.linalg.norm(x_temp - point) < tol:
                return False
        
        # Проверяем, что период срабатывает
        x_period = self.period_mapping(point)
        return np.linalg.norm(x_period - point) < tol


# Пример использования с отображением Хенона в 3D
# def LermanMap(x: np.ndarray, eps: float = 0.0505, alpha: float = 0.431, beta: float = 0.2, p: float = 1) -> np.ndarray:
#     """
#     3D отображение Lerman.
    
#     Parameters:
#     -----------
#     x : np.ndarray
#         Точка в R^3 [ksi, eta, teta]
#     eps, alpha, beta, p : float
#         Параметры отображения
        
#     Returns:
#     --------
#     next_x : np.ndarray
#         Следующая точка
#     """
#     # LOCAL 2
#     ro = np.sqrt(x[0]**2+x[1]**2)/p
#     x[2] = x[2]
#     phi0 = np.atan2(x[1], x[0]) + x[2]

#     r = p * (ro / p) ** ((alpha+eps)/(alpha-eps))
#     teta1 = ((beta + eps) / (alpha - eps)) * np.log(p/ro) + x[2]
#     phi1 = ((beta - eps) / (alpha - eps)) * np.log(p/ro) + phi0
#     # phi1 = ((((beta - eps) / (alpha - eps)) * np.log(p/ro) + phi0) % (2 * np.pi)) - np.pi
#     while (phi1 > np.pi):
#         phi1 = phi1 - 2 * np.pi
#     while (phi1 < -np.pi):
#         phi1 += 2 * np.pi

#     ksi1 = r * p * np.cos(phi1-teta1)
#     eta1 = r * p * np.sin(phi1-teta1)
#     phi1 = phi1

#     # ******** S3
#     ksi2 = ksi1 + eps * np.cos(2*phi1)
#     eta2 = eta1 + eps * np.sin(2*phi1)
#     # teta2 = phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
#     teta2 = phi1
#     while teta2>np.pi:
#         teta2 -= 2 * np.pi
#     while teta2<-np.pi:
#         teta2 += 2 * np.pi

#     return np.array([ksi2, eta2, teta2])

def LermanMap(x: np.ndarray, eps: float = 0.0505, alpha: float = 0.431, beta: float = 0.2, p: float = 1):
    # eps, alpha, beta, p = params
    # print(eps, alpha, p_s, p_u)
    ksi0, eta0, teta0 = x

    # LOCAL 1
    # r_factor = (p**2) / np.sqrt(ksi0**2 + eta0**2)
    # scale = (r_factor) ** ((-2 * eps) / (alpha - eps))

    # Phi0 = np.arctan2(eta0, ksi0)
    # # Phi1 = (-(2 * eps) / (alpha - eps)) * np.log(r_factor) + Phi0
    # # if Phi1>np.pi:
    # #     Phi1 -= 2 * np.pi
    # # elif Phi1<-np.pi:
    # #     Phi1 += 2 * np.pi

    # # ksi1 = p_s * p_u * scale * np.cos(Phi1)
    # # eta1 = p_s * p_u * scale * np.sin(Phi1)
    # # phi1 = Phi1 - teta0
    # ksi1 = scale * ( ksi0 * np.cos( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.sin(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # eta1 = scale * (-ksi0 * np.sin( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.cos(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # # phi1 = (teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)) % (2*np.pi)
    # phi1 =  teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)

    # LOCAL 2
    ro = np.sqrt(ksi0**2+eta0**2)/p
    teta0 = teta0
    phi0 = np.atan2(eta0, ksi0) + teta0

    # while (phi0 > m.pi):
    #     phi0 -= 2 * m.pi
    # while (phi0 < -m.pi):
    #     phi0 += 2 * m.pi

    # print("ro, teta0, phi0")
    # print(ro, teta0, phi0)

    r = p * (ro / p) ** ((alpha+eps)/(alpha-eps))
    teta1 = ((beta + eps) / (alpha - eps)) * np.log(p/ro) + teta0
    phi1 = (((beta - eps) / (alpha - eps)) * np.log(p/ro) + phi0)

    while (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    while (phi1 < -m.pi):
        phi1 += 2 * m.pi
    # print("r, teta1, phi1")
    # print(r, teta1, phi1)
    ksi1 = r * p * np.cos(phi1-teta1)
    eta1 = r * p * np.sin(phi1-teta1)
    phi1 = phi1

    # LOCAL 3
    # r_factor = (p * p)
    # scale = pow(r_factor, (-2.0 * eps) / (alpha - eps))
    # Phi0 = m.atan2(eta0, ksi0)
    # # angle = (2.0 * eps / (alpha - eps)) * m.log(r_factor)
        
    # ksi1 = scale * ksi0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # eta1 = scale * eta0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # phi1 = teta0 + Phi0 + ((beta - eps) / (alpha - eps)) * m.log(r_factor/(np.sqrt(ksi0**2+eta0**2)))
    # while (phi1 > m.pi):
    #     phi1 = phi1 - 2 * m.pi
    # while (phi1 < -m.pi):
    #     phi1 += 2 * m.pi

    # ******** S3
    # ksi2 = ksi1 + eps * np.cos(2*phi1)
    # eta2 = eta1 + eps * np.sin(2*phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S2
    # ksi2 = ksi1 + eps * np.cos(phi1)
    # eta2 = eta1 + eps * np.sin(phi1)
    # # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S1
    ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    teta2 =  phi1
    while teta2>np.pi:
        teta2 -= 2 * np.pi
    while teta2<-np.pi:
        teta2 += 2 * np.pi

    return np.array([ksi2,eta2,teta2])


def example_3d_henon(x: np.ndarray, a: float = 1.4, b: float = 0.3, c: float = 0.1) -> np.ndarray:
    x_new = 1 - a * x[0]**2 + x[1]
    y_new = b * x[0]
    z_new = c * x[1] + (1 - c) * x[2]
    return np.array([x_new, y_new, z_new])

# fig, ax = plt.subplots()
# x = np.array([0.001,0.001,0.001])
# points = []
# for i in range(155000):
#     x = LermanMap(x)
#     if i>105000:
#         points.append(x.copy())
# points = np.array(points).T
# print(points.shape)
# print(points)
# for i in range(100):
#     print(points[2][i])
# ax.plot(points[2][::10], points[1][::10], marker="o", markersize=0.3, linestyle="", color="black")
# plt.show()

if __name__ == "__main__":
    print("Пример 1: Поиск периодической точки для LermanMap")
    print("-" * 60)
    
    # Создаем отображение с фиксированными параметрами
    per = 10
    def map(x):
        return LermanMap(x, eps=0.014, alpha=0.048, beta=0.2, p=1)
    
    # Инициализируем поисковики под разные периоды
    finders = []
    for i in range(per):
        finder = PeriodicPointFinder3D(map, period=i+1)
        finders.append(finder)
    
    # Начальные приближения ksi eta teta
    initial_guesses = [
        #[0.05, 0.43, 0.2, 1]
        # np.array([-0.004, -0.001067, -1.1867])
        # np.array([-0.00374156, -0.0009989, 2.053225])
        # np.array([0.08049354, -0.04534888, 0.6200527])
        # np.array([-0.08448109, 0.03617788, -1.005346])
        np.array([0.009804907368803155, 0.005050402256139094, 0.49902186817420713], dtype=np.longdouble)
        # np.array([-0.0010371836, 0.0067190791, 1.9652417036])
        # 0.071253 -0.002914 -0.006253
        #[0.0505, 0.431, 0.2, 1]
        # np.array([-0.00577228, -0.00222498, 2.0892851])
        # np.array([0.08212919, -0.04518525, -2.4712331]),
        # np.array([0.0341986412, 0.0177492215, -1.9813259076])
    ]
    
    # Ищем периодические точки
    results = []
    for i in range(per):
        result = finders[i].find_multiple_points(initial_guesses, tol=1e-12)
        # print(result)
        results.append(result[0])

    # Выводим результат
    for i, res in enumerate(results):
        # res = result[0]
        print(f"\nРезультат {i+1} (начальное приближение: {initial_guesses}):")
        # print(f"\nРезультат {i+1} (начальное приближение: {initial_guesses[0]:20f}, {initial_guesses[1]:20f}, {initial_guesses[2]:20f}):")
        if res['success']:
            print(f"  Найдена точка: {res['point']}")
            print(f"  Расстояние от приближения: {res['distance']}")
            print(f"  Остаток: {res['residual']:.2e}")
            print(f"  Собственные числа: {res['eigenvalues']}")
            print(f"  Тип устойчивости: {res['stability']}")
            print(f"  Сообщение: {res['message']}")
            print(f"  iterations: {res['iterations']}")
        else:
            print(f"  Не удалось найти точку: {res['message']}")
    
    # # Пример 3: Анализ производной в точке
    # print("\n" + "=" * 60)
    # print("Пример 3: Проверка численного якобиана")
    # print("-" * 60)
    
    # # Тестовая точка
    # test_point = np.array([0.5, 0.2, 0.3])
    
    # # Якобиан в точке
    # J = finder.numerical_jacobian(test_point)
    # print(f"Якобиан в точке {test_point}:")
    # print(J)
    
    # # Проверим точность дифференцирования
    # print(f"\nПроверка точности (eps = {finder.eps}):")
    
    # # Аналитический якобиан для 3D Хенона (если известен)
    # # Для общего случая проверим выполнение условия симметрии
    # J_T = finder.numerical_jacobian(test_point + 1e-6)
    # diff = np.linalg.norm(J_T - J)
    # print(f"Изменение якобиана при малом сдвиге: {diff:.2e}")