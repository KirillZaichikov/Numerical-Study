import sys
import os

# Добавляем родительскую директорию в sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Теперь можно импортировать модуль из родительской директории
from utils.integrator import *

import numpy as np
import matplotlib.pyplot as plt
from typing import Callable, Tuple, Optional, List
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
import time

# Игнорируем предупреждения о делении на ноль в методе Ньютона
warnings.filterwarnings('ignore')
np.set_printoptions(suppress=True)

class NewtonSolver:
    """Класс для решения нелинейных уравнений методом Ньютона"""
    
    def __init__(self, dim: int = 5, tol: float = 0.0001, max_iter: int = 50, eps: float = 1e-8):
        """
        Инициализация решателя Ньютона
        
        Parameters:
        -----------
        dim : int
            Размерность пространства (по умолчанию 5)
        tol : float
            Допуск для сходимости
        max_iter : int
            Максимальное число итераций
        eps : float
            Малое число для вычисления конечно-разностного якобиана
        """
        self.dim = dim
        self.tol = tol
        self.max_iter = max_iter
        self.eps = eps
    
    def finite_difference_jacobian(self, G: Callable, x: np.ndarray, 
                                 alpha: float, beta: float) -> np.ndarray:
        """
        Вычисление якобиана методом конечных разностей
        
        Parameters:
        -----------
        G : Callable
            Функция G(x, alpha, beta) = F(x, alpha, beta) - x
        x : np.ndarray
            Точка, в которой вычисляется якобиан
        alpha : float
            Значение параметра alpha
        beta : float
            Значение параметра beta
            
        Returns:
        --------
        J : np.ndarray
            Матрица Якоби размерности (dim, dim)
        """
        J = np.zeros((self.dim, self.dim))
        G0 = G(x, alpha, beta)
        
        for j in range(self.dim):
            # Создаем возмущенный вектор
            dx = np.zeros(self.dim)
            dx[j] = self.eps
            # pert = fix_integrals(alpha, beta, x + dx)
            # print(pert)
            G_perturbed = G(x+dx, alpha, beta)
            # print("G0",G0)
            # print("G_pert",G_perturbed)
            J[:, j] = (G_perturbed - G0) / self.eps
            
        return J
    
    def solve(self, G: Callable, x0: np.ndarray, 
              alpha: float, beta: float) -> Tuple[Optional[np.ndarray], int, bool]:
        """
        Решение уравнения G(x, alpha, beta) = 0 методом Ньютона
        
        Parameters:
        -----------
        G : Callable
            Функция G(x, alpha, beta) = F(x, alpha, beta) - x
        x0 : np.ndarray
            Начальное приближение
        alpha : float
            Значение параметра alpha
        beta : float
            Значение параметра beta
            
        Returns:
        --------
        x_sol : Optional[np.ndarray]
            Найденное решение или None, если не сошлось
        iterations : int
            Число выполненных итераций
        success : bool
            Флаг успешности сходимости
        """
        x = x0.copy()
        
        for iteration in range(self.max_iter):
            # Вычисляем G(x)
            # print(f"Номер итерации {iteration}")
            G_val = G(x, alpha, beta)
            
            # Проверяем сходимость
            # print(np.linalg.norm(G_val))
            if np.linalg.norm(G_val) < self.tol:
                return x, iteration + 1, True
            
            # Вычисляем якобиан
            # print(f"Считаем Якобиан {alpha}")
            J = self.finite_difference_jacobian(G, x, alpha, beta)
            # print(f"Посчитали Якобиан {alpha}")
            
            # Решаем линейную систему J * dx = -G(x)
            try:
                # print(f"Решаем лин систему {alpha}")
                dx = np.linalg.solve(J, -G_val)
                # print(dx)
                # print(f"Решили лин систему {alpha}")
            except np.linalg.LinAlgError:
                # Если матрица вырожденная, пытаемся использовать псевдообратную
                try:
                    dx = -np.linalg.pinv(J) @ G_val
                except:
                    return None, iteration + 1, False
            
            # Обновляем решение
            x += dx
            # x = fix_integrals(alpha, beta, x)
            # print(x)
            # print(x)
            # Проверяем на расходимость (очень большие значения)
            if np.any(np.abs(x) > 1e10):
                return None, iteration + 1, False
        
        # Если вышли за пределы максимального числа итераций
        G_val = G(x, alpha, beta)
        if np.linalg.norm(G_val) < self.tol:
            return x, self.max_iter, True
        else:
            return None, self.max_iter, False


def process_point_worker(args):
    obj, point, beta_min, beta_max, beta_step = args

    print("PID:", os.getpid())

    alpha = point['alpha']
    x_alpha = point['x']
    beta0 = point['beta']

    down_branch, up_branch, down_branch_MG, up_branch_MG = obj.vertical_continuation_for_alpha(
        x_alpha, alpha, beta0, beta_min, beta_max, beta_step
    )

    full_line = []
    full_line_MG = []

    if down_branch:
        down_branch_sorted = sorted(down_branch, key=lambda p: p['beta'], reverse=True)
        full_line.extend(down_branch_sorted)
        down_branch_MG_sorted = sorted(down_branch_MG, key=lambda p: p['beta'], reverse=True)
        full_line_MG.extend(down_branch_MG_sorted)

    if point['beta'] == beta0 and (not full_line or full_line[-1]['beta'] != beta0):
        full_line.append(point)
        full_line_MG.append(point)

    if up_branch:
        up_branch_sorted = sorted(up_branch, key=lambda p: p['beta'])
        full_line.extend(up_branch_sorted)
        up_branch_MG_sorted = sorted(up_branch_MG, key=lambda p: p['beta'])
        full_line_MG.extend(up_branch_MG_sorted)

    return full_line, full_line_MG


class FixedPointContinuation:
    """Класс для протягивания неподвижных точек по параметрам"""
    
    def __init__(self, F: Callable, dim: int = 5):
        """
        Инициализация класса
        
        Parameters:
        -----------
        F : Callable
            Функция F(x, alpha, beta)
        dim : int
            Размерность пространства (по умолчанию 5)
        """
        self.F = F
        self.dim = dim
        self.solver = NewtonSolver(dim=dim)
        
    def G(self, x: np.ndarray, alpha: float, beta: float) -> np.ndarray:
        """
        Функция G(x, alpha, beta) = F(x, alpha, beta) - x
        
        Parameters:
        -----------
        x : np.ndarray
            Вектор состояния
        alpha : float
            Первый параметр
        beta : float
            Второй параметр
            
        Returns:
        --------
        np.ndarray
            Значение функции G
        """
        return self.F(x, alpha, beta) - x
    
    def horizontal_continuation(self, x_start: np.ndarray, 
                               alpha0: float, beta0: float,
                               alpha_min: float, alpha_max: float,
                               alpha_step: float) -> Tuple[List[dict], List[dict]]:
        """
        Горизонтальное протягивание вдоль alpha при фиксированном beta
        
        Parameters:
        -----------
        x_start : np.ndarray
            Начальное приближение при (alpha0, beta0)
        alpha0 : float
            Начальное значение alpha
        beta0 : float
            Фиксированное значение beta
        alpha_min : float
            Минимальное значение alpha
        alpha_max : float
            Максимальное значение alpha
        alpha_step : float
            Шаг по alpha
            
        Returns:
        --------
        left_branch : List[dict]
            Ветвь влево (alpha уменьшается)
        right_branch : List[dict]
            Ветвь вправо (alpha увеличивается)
        """
        # Ветвь влево (alpha уменьшается)
        left_branch = []
        left_branch_MG = []
        current_x = x_start.copy()
        current_alpha = alpha0
        
        # Сохраняем начальную точку
        # print("поиск первой точки")
        sol, iterations, success = self.solver.solve(self.G, current_x, alpha0, beta0)
        if success:
            current_x = sol
            left_branch.append({
                'alpha': alpha0,
                'beta': beta0,
                'x': sol.copy(),
                'iterations': iterations,
                'success': success
            })

            left_branch_MG.append({
                'alpha': alpha0,
                'beta': beta0,
                'x': ADtoMG(sol.copy(), alpha0, beta0).copy(),
                'iterations': iterations,
                'success': success
            })
        # print("нашли", sol, "за", iterations, "итераций")

        # Протягиваем влево
        alpha = alpha0 - alpha_step
        while alpha >= alpha_min - 1e-10:
            # Используем решение с предыдущего шага как начальное приближение
            print(f"Поиск точки {alpha}")
            # fix_integrals(alpha, beta0, current_x)
            sol, iterations, success = self.solver.solve(self.G, current_x, alpha, beta0)
            # print("нашли", sol, "за", iterations, "итераций")
            
            if success:
                current_x = sol
                left_branch.append({
                    'alpha': alpha,
                    'beta': beta0,
                    'x': sol.copy(),
                    'iterations': iterations,
                    'success': success
                })
                left_branch_MG.append({
                    'alpha': alpha,
                    'beta': beta0,
                    'x': ADtoMG(sol.copy(), alpha, beta0).copy(),
                    'iterations': iterations,
                    'success': success
                })
                alpha -= alpha_step
            else:
                # Если не сошлось, прерываем продолжение
                break
        
        # Ветвь вправо (alpha увеличивается)
        right_branch = []
        right_branch_MG = []
        current_x = x_start.copy()
        current_alpha = alpha0
        
        # Уже добавили начальную точку, но нужно проверить сходимость
        if len(left_branch) == 0:
            # Если начальная точка не сошлась, пытаемся снова
            sol, iterations, success = self.solver.solve(self.G, current_x, alpha0, beta0)
            if success:
                current_x = sol
                right_branch.append({
                    'alpha': alpha0,
                    'beta': beta0,
                    'x': sol.copy(),
                    'iterations': iterations,
                    'success': success
                })
                right_branch_MG.append({
                    'alpha': alpha0,
                    'beta': beta0,
                    'x': ADtoMG(sol.copy(), alpha0, beta0).copy(),
                    'iterations': iterations,
                    'success': success
                })
        
        # Протягиваем вправо
        alpha = alpha0 + alpha_step
        while alpha <= alpha_max + 1e-10:
            sol, iterations, success = self.solver.solve(self.G, current_x, alpha, beta0)
            
            if success:
                current_x = sol
                right_branch.append({
                    'alpha': alpha,
                    'beta': beta0,
                    'x': sol.copy(),
                    'iterations': iterations,
                    'success': success
                })
                right_branch_MG.append({
                    'alpha': alpha,
                    'beta': beta0,
                    'x': ADtoMG(sol.copy(), alpha, beta0).copy(),
                    'iterations': iterations,
                    'success': success
                })
                alpha += alpha_step
            else:
                break
        
        return left_branch, right_branch, left_branch_MG, right_branch_MG
    
    def vertical_continuation_for_alpha(self, x_alpha: np.ndarray, 
                                       alpha: float, beta0: float,
                                       beta_min: float, beta_max: float,
                                       beta_step: float) -> Tuple[List[dict], List[dict]]:
        """
        Вертикальное протягивание вдоль beta для фиксированного alpha
        
        Parameters:
        -----------
        x_alpha : np.ndarray
            Начальное приближение при (alpha, beta0)
        alpha : float
            Фиксированное значение alpha
        beta0 : float
            Начальное значение beta
        beta_min : float
            Минимальное значение beta
        beta_max : float
            Максимальное значение beta
        beta_step : float
            Шаг по beta
            
        Returns:
        --------
        down_branch : List[dict]
            Ветвь вниз (beta уменьшается)
        up_branch : List[dict]
            Ветвь вверх (beta увеличивается)
        """
        # Ветвь вниз (beta уменьшается)
        down_branch = []
        down_branch_MG = []
        current_x = x_alpha.copy()
        
        # Проверяем начальную точку
        sol, iterations, success = self.solver.solve(self.G, current_x, alpha, beta0)
        if success:
            current_x = sol
            down_branch.append({
                'alpha': alpha,
                'beta': beta0,
                'x': sol.copy(),
                'iterations': iterations,
                'success': success
            })
            down_branch_MG.append({
                'alpha': alpha,
                'beta': beta0,
                'x': ADtoMG(sol.copy(), alpha, beta0).copy(),
                'iterations': iterations,
                'success': success
            })
        
        # Протягиваем вниз
        beta = beta0 - beta_step
        while beta >= beta_min - 1e-10:
            sol, iterations, success = self.solver.solve(self.G, current_x, alpha, beta)
            
            if success:
                current_x = sol
                down_branch.append({
                    'alpha': alpha,
                    'beta': beta,
                    'x': sol.copy(),
                    'iterations': iterations,
                    'success': success
                })
                down_branch_MG.append({
                    'alpha': alpha,
                    'beta': beta,
                    'x': ADtoMG(sol.copy(), alpha, beta).copy(),
                    'iterations': iterations,
                    'success': success
                })
                beta -= beta_step
            else:
                break
        
        # Ветвь вверх (beta увеличивается)
        up_branch = []
        up_branch_MG = []
        current_x = x_alpha.copy()
        
        if len(down_branch) == 0:
            # Если начальная точка не сошлась, пытаемся снова
            sol, iterations, success = self.solver.solve(self.G, current_x, alpha, beta0)
            if success:
                current_x = sol
                up_branch.append({
                    'alpha': alpha,
                    'beta': beta0,
                    'x': sol.copy(),
                    'iterations': iterations,
                    'success': success
                })
                up_branch_MG.append({
                    'alpha': alpha,
                    'beta': beta0,
                    'x': ADtoMG(sol.copy(), alpha, beta0).copy(),
                    'iterations': iterations,
                    'success': success
                })
        
        # Протягиваем вверх
        beta = beta0 + beta_step
        while beta <= beta_max + 1e-10:
            sol, iterations, success = self.solver.solve(self.G, current_x, alpha, beta)
            
            if success:
                current_x = sol
                up_branch.append({
                    'alpha': alpha,
                    'beta': beta,
                    'x': sol.copy(),
                    'iterations': iterations,
                    'success': success
                })
                up_branch_MG.append({
                    'alpha': alpha,
                    'beta': beta,
                    'x': ADtoMG(sol.copy(), alpha, beta).copy(),
                    'iterations': iterations,
                    'success': success
                })
                beta += beta_step
            else:
                break
        
        return down_branch, up_branch, down_branch_MG, up_branch_MG
    
    """
    # def parallel_vertical_continuation(self, horizontal_points: List[dict],
    #                                   beta_min: float, beta_max: float,
    #                                   beta_step: float, n_workers: int = 12) -> List[List[dict]]:
    #     
    #     Параллельное вертикальное протягивание для всех точек горизонтальной линии
        
    #     Parameters:
    #     -----------
    #     horizontal_points : List[dict]
    #         Точки горизонтальной линии
    #     beta_min : float
    #         Минимальное значение beta
    #     beta_max : float
    #         Максимальное значение beta
    #     beta_step : float
    #         Шаг по beta
    #     n_workers : int
    #         Число потоков для параллельных вычислений
            
    #     Returns:
    #     --------
    #     vertical_lines : List[List[dict]]
    #         Список вертикальных линий для каждого alpha
    #     
    #     vertical_lines = []
        
    #     # Функция для обработки одной точки горизонтальной линии
    #     def process_point(point):
    #         alpha = point['alpha']
    #         x_alpha = point['x']
    #         beta0 = point['beta']
            
    #         down_branch, up_branch = self.vertical_continuation_for_alpha(
    #             x_alpha, alpha, beta0, beta_min, beta_max, beta_step
    #         )
            
    #         # Объединяем ветви (вниз + начальная + вверх)
    #         # Убираем дубликаты начальной точки
    #         full_line = []
    #         if down_branch:
    #             # Сортируем по убыванию beta
    #             down_branch_sorted = sorted(down_branch, key=lambda p: p['beta'], reverse=True)
    #             full_line.extend(down_branch_sorted)
            
    #         # Добавляем начальную точку, если ее еще нет
    #         if point['beta'] == beta0 and (not full_line or full_line[-1]['beta'] != beta0):
    #             full_line.append(point)
            
    #         if up_branch:
    #             # Сортируем по возрастанию beta
    #             up_branch_sorted = sorted(up_branch, key=lambda p: p['beta'])
    #             full_line.extend(up_branch_sorted)
            
    #         return full_line
        
    #     # Параллельная обработка
    #     # with ThreadPoolExecutor(max_workers=n_workers) as executor:
    #     #     futures = [executor.submit(process_point, point) for point in horizontal_points]
            
    #     #     for future in as_completed(futures):
    #     #         vertical_lines.append(future.result())

    #     with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
    #         futures = [executor.submit(process_point, point) 
    #            for point in horizontal_points]

    #         for future in as_completed(futures):
    #             vertical_lines.append(future.result())
        
    #     return vertical_lines
    """

    def parallel_vertical_continuation(self, horizontal_points, beta_min, beta_max, beta_step, n_workers=10):
        args = [
            (self, point, beta_min, beta_max, beta_step)
            for point in horizontal_points
        ]

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            vertical_lines = list(executor.map(process_point_worker, args))
            vl = []
            vl_MG = []
            for elements in vertical_lines:
                vl.append(elements[0])
                vl_MG.append(elements[1])


        return vl, vl_MG

    def full_continuation(self, x_start: np.ndarray, 
                         alpha0: float, beta0: float,
                         alpha_min: float, alpha_max: float, alpha_step: float,
                         beta_min: float, beta_max: float, beta_step: float,
                         n_workers: int = 8) -> dict:
        """
        Полное протягивание по двум параметрам
        
        Parameters:
        -----------
        x_start : np.ndarray
            Начальное приближение
        alpha0, beta0 : float
            Начальные значения параметров
        alpha_min, alpha_max, alpha_step : float
            Диапазон и шаг по alpha
        beta_min, beta_max, beta_step : float
            Диапазон и шаг по beta
        n_workers : int
            Число потоков для параллельных вычислений
            
        Returns:
        --------
        result : dict
            Словарь с результатами
        """
        print("Фаза 1: Горизонтальное протягивание...")
        start = time.time()
        left_branch, right_branch, left_branch_MG, right_branch_MG = self.horizontal_continuation(
            x_start, alpha0, beta0, alpha_min, alpha_max, alpha_step
        )
        print("Время на горизонтальное протягивание", time.time()-start)
        
        # Объединяем горизонтальную линию
        horizontal_line = []
        if left_branch:
            # Сортируем по возрастанию alpha
            left_branch_sorted = sorted(left_branch, key=lambda p: p['alpha'])
            horizontal_line.extend(left_branch_sorted)
        
        # Добавляем правую ветвь
        if right_branch:
            horizontal_line.extend(right_branch)
        
        print(f"Найдено {len(horizontal_line)} точек на горизонтальной линии")
        
        print("Фаза 2: Вертикальное протягивание...")
        start = time.time()
        vertical_lines, vertical_lines_MG = self.parallel_vertical_continuation(
            horizontal_line, beta_min, beta_max, beta_step, n_workers
        )
        print("Время на вертикальное протягивание", time.time()-start)
        
        # Собираем все точки в сетку
        all_points = []
        for line in vertical_lines:
            all_points.extend(line)
        all_points_MG = []
        for line in vertical_lines_MG:
            all_points_MG.extend(line)
        print(f"Всего найдено {len(all_points)} точек c неподвижной точкой")
        print(f"Всего найдено {len(all_points_MG)} в M gamma")
        
        return {
            'horizontal_line': horizontal_line,
            'vertical_lines': vertical_lines,
            'all_points': all_points,
            'all_points_MG': all_points_MG
        }

'''
# Тестовая функция F
# def test_F(x: np.ndarray, alpha: float, beta: float) -> np.ndarray:
#     """
#     Тестовая функция F(x, alpha, beta)
#     """
#     B=0.7
#     dim = len(x)
#     M_1 = alpha
#     M_2 = beta
#     result = np.zeros(dim)
#     result[0] = x[1]
#     result[1] = x[2]
#     result[2] = M_1 + B * x[0] + M_2 * x[1] - x[2]**2
#     # for i in range(dim):
#     #     # Циклическое условие
#     #     next_idx = (i + 1) % dim
#     #     result[i] = (1 - alpha) * x[i] + beta * np.sin(x[next_idx])
    
#     return result
'''

def calc_r_vec(a1, a2, h, gamma, d):
    gamma1, gamma2, gamma3 = gamma[0], gamma[1], gamma[2]
    r1 = -(a1 * gamma1) / gamma3 # КУЗНЕЦОВ
    r2 = -(a2 * gamma2) / gamma3 # КУЗНЕЦОВ
    r3 = -h + (1/2)*((a1 * gamma1 ** 2 + a2 * gamma2 ** 2)/(gamma3 ** 2)) # КУЗНЕЦОВ
    r = np.array([r1,r2,r3])
    return r

def calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M):
    Q = np.array([[np.cos(d), np.sin(d), 0],
                 [-np.sin(d),np.cos(d),0],
                 [0,0,1]])
    Q = np.array([[np.cos(d), np.sin(d), 0],
                 [-np.sin(d),np.cos(d),0],
                 [0,0,1]])
    # Формула бизяева в maple
    A = Q @ np.diag([I1,I2,I3]) @ Q.T + np.identity(3) * np.dot(r,r) - np.outer(r, r)
    A_inv = np.linalg.inv(A)
    # print("JQ_inv\n", A_inv)
    omega = A_inv @ M
    return omega

# Нужны для следующих двух функций
params = np.array([0.485, 2, 6, 7, 9, 4, 1, 752, 100])

def ADtoMG(x, alpha, beta):
    params[0] = alpha
    params[-2] = beta

    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    ### Замена от переменных Андуайе-Депри
    l = x[0]
    LG = x[1]
    HG = x[2]

    m1 = (1-LG**2)**(1/2) * np.sin(l)
    m2 = (1-LG**2)**(1/2) * np.cos(l)
    m3 = LG
    M = np.array([m1,m2,m3]) # При такой замене вектор M будет единичной длины

    gamma1 = (HG * (1-LG**2)**(1/2) + LG*(1-HG**2)**(1/2)) * np.sin(l)
    gamma2 = (HG * (1-LG**2)**(1/2) + LG*(1-HG**2)**(1/2)) * np.cos(l)
    gamma3 = HG * LG - (1-LG**2)**(1/2) * (1-HG**2)**(1/2) * 1
    gamma = np.array([gamma1,gamma2,gamma3]) # гамма тоже единичный
    gamma = -gamma # надо повернуть

    QC = np.array([[np.cos(d), np.sin(d), 0],
                [-np.sin(d),np.cos(d),0],
                [0,0,1]])
    M = QC @ M
    gamma = QC @ gamma

    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    tmp = M.copy()
    for k in range(3):
        M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))
    fx = np.array([*M, *gamma])
    return fx

def diffFunc(x, fx, params, H):
    delta = params[0]
    J = np.array(params[1:4])
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    M = np.array(x[:3])
    gamma = np.array(x[3:])
    # print(M)
    # print(gamma)

    r = np.ndarray(3)
    r[0] = -params[4] * gamma[0] / gamma[2]
    r[1] = -params[5] * gamma[1] / gamma[2]
    r[2] = -params[6] + (params[4] * pow(gamma[0], 0.2e1) + params[5] * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.2e1) / 0.2e1
    
    JQ_kuz = np.ndarray((3,3))

    JQ_kuz[0][0] = pow(np.cos(delta), 0.2e1) * J[0] + pow(np.sin(delta), 0.2e1) * J[1] + pow(r[1], 0.2e1) + pow(r[2], 0.2e1)
    JQ_kuz[0][1] = -np.cos(delta) * J[0] * np.sin(delta) + np.sin(delta) * J[1] * np.cos(delta) - r[0] * r[1]
    JQ_kuz[0][2] = -r[0] * r[2]
    JQ_kuz[1][0] = -np.cos(delta) * J[0] * np.sin(delta) + np.sin(delta) * J[1] * np.cos(delta) - r[0] * r[1]
    JQ_kuz[1][1] = pow(np.sin(delta), 0.2e1) * J[0] + pow(np.cos(delta), 0.2e1) * J[1] + pow(r[0], 0.2e1) + pow(r[2], 0.2e1)
    JQ_kuz[1][2] = -r[1] * r[2]
    JQ_kuz[2][0] = -r[0] * r[2]
    JQ_kuz[2][1] = -r[1] * r[2]
    JQ_kuz[2][2] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) + J[2]

    JQ_rev = np.linalg.inv(JQ_kuz)

    omega = JQ_rev @ M

    dgamma = np.ndarray(3)
    dgamma[0] =  gamma[1] * omega[2] - gamma[2] * omega[1]
    dgamma[1] = -gamma[0] * omega[2] + gamma[2] * omega[0]
    dgamma[2] =  gamma[0] * omega[1] - gamma[1] * omega[0]

    Jr = np.ndarray((3,3))

    Jr[0][0] = -a1 / gamma[2]
    Jr[0][1] = 0
    Jr[0][2] = a1 * gamma[0] * pow(gamma[2], -0.2e1)
    Jr[1][0] = 0
    Jr[1][1] = -a2 / gamma[2]
    Jr[1][2] = a2 * gamma[1] * pow(gamma[2], -0.2e1)
    Jr[2][0] = a1 * gamma[0] * pow(gamma[2], -0.2e1)
    Jr[2][1] = a2 * gamma[1] * pow(gamma[2], -0.2e1)
    Jr[2][2] = -(a1 * pow(gamma[0], 0.2e1) + a2 * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.3e1)

    dr = np.ndarray(3)

    dr[0] = Jr[0][0] * dgamma[0] + Jr[0][1] * dgamma[1] + Jr[0][2] * dgamma[2]
    dr[1] = Jr[1][0] * dgamma[0] + Jr[1][1] * dgamma[1] + Jr[1][2] * dgamma[2]
    dr[2] = Jr[2][0] * dgamma[0] + Jr[2][1] * dgamma[1] + Jr[2][2] * dgamma[2]

    dM = np.ndarray(3)

    dM[0] =  M[1] * omega[2] - M[2] * omega[1] + dr[1] * ( omega[0] * r[1] - omega[1] * r[0]) - dr[2] * (-omega[0] * r[2] + omega[2] * r[0]) + g0 * (-r[2] * gamma[1] + r[1] * gamma[2])
    dM[1] = -M[0] * omega[2] + M[2] * omega[0] - dr[0] * ( omega[0] * r[1] - omega[1] * r[0]) + dr[2] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * ( r[2] * gamma[0] - r[0] * gamma[2])
    dM[2] =  M[0] * omega[1] - M[1] * omega[0] + dr[0] * (-omega[0] * r[2] + omega[2] * r[0]) - dr[1] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * (-r[1] * gamma[0] + r[0] * gamma[1])

    fx[0] = dM[0]
    fx[1] = dM[1]
    fx[2] = dM[2]
    fx[3] = dgamma[0]
    fx[4] = dgamma[1]
    fx[5] = dgamma[2]

def diffFuncP(x, fx, params, H):
    delta = params[0]
    J = np.array(params[1:4])
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    M = np.array(x[:3])
    gamma = np.array(x[3:])

    r = np.ndarray(3)
    r[0] = -params[4] * gamma[0] / gamma[2]
    r[1] = -params[5] * gamma[1] / gamma[2]
    r[2] = -params[6] + (params[4] * pow(gamma[0], 0.2e1) + params[5] * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.2e1) / 0.2e1
    
    JQ_kuz = np.ndarray((3,3))

    JQ_kuz[0][0] = pow(np.cos(delta), 0.2e1) * J[0] + pow(np.sin(delta), 0.2e1) * J[1] + pow(r[1], 0.2e1) + pow(r[2], 0.2e1)
    JQ_kuz[0][1] = -np.cos(delta) * J[0] * np.sin(delta) + np.sin(delta) * J[1] * np.cos(delta) - r[0] * r[1]
    JQ_kuz[0][2] = -r[0] * r[2]
    JQ_kuz[1][0] = -np.cos(delta) * J[0] * np.sin(delta) + np.sin(delta) * J[1] * np.cos(delta) - r[0] * r[1]
    JQ_kuz[1][1] = pow(np.sin(delta), 0.2e1) * J[0] + pow(np.cos(delta), 0.2e1) * J[1] + pow(r[0], 0.2e1) + pow(r[2], 0.2e1)
    JQ_kuz[1][2] = -r[1] * r[2]
    JQ_kuz[2][0] = -r[0] * r[2]
    JQ_kuz[2][1] = -r[1] * r[2]
    JQ_kuz[2][2] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) + J[2]

    JQ_rev = np.linalg.inv(JQ_kuz)

    omega = JQ_rev @ M

    dgamma = np.ndarray(3)
    dgamma[0] =  gamma[1] * omega[2] - gamma[2] * omega[1]
    dgamma[1] = -gamma[0] * omega[2] + gamma[2] * omega[0]
    dgamma[2] =  gamma[0] * omega[1] - gamma[1] * omega[0]

    Jr = np.ndarray((3,3))

    Jr[0][0] = -a1 / gamma[2]
    Jr[0][1] = 0
    Jr[0][2] = a1 * gamma[0] * pow(gamma[2], -0.2e1)
    Jr[1][0] = 0
    Jr[1][1] = -a2 / gamma[2]
    Jr[1][2] = a2 * gamma[1] * pow(gamma[2], -0.2e1)
    Jr[2][0] = a1 * gamma[0] * pow(gamma[2], -0.2e1)
    Jr[2][1] = a2 * gamma[1] * pow(gamma[2], -0.2e1)
    Jr[2][2] = -(a1 * pow(gamma[0], 0.2e1) + a2 * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.3e1)

    dr = np.ndarray(3)

    dr[0] = Jr[0][0] * dgamma[0] + Jr[0][1] * dgamma[1] + Jr[0][2] * dgamma[2]
    dr[1] = Jr[1][0] * dgamma[0] + Jr[1][1] * dgamma[1] + Jr[1][2] * dgamma[2]
    dr[2] = Jr[2][0] * dgamma[0] + Jr[2][1] * dgamma[1] + Jr[2][2] * dgamma[2]

    dM = np.ndarray(3)

    dM[0] =  M[1] * omega[2] - M[2] * omega[1] + dr[1] * ( omega[0] * r[1] - omega[1] * r[0]) - dr[2] * (-omega[0] * r[2] + omega[2] * r[0]) + g0 * (-r[2] * gamma[1] + r[1] * gamma[2])
    dM[1] = -M[0] * omega[2] + M[2] * omega[0] - dr[0] * ( omega[0] * r[1] - omega[1] * r[0]) + dr[2] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * ( r[2] * gamma[0] - r[0] * gamma[2])
    dM[2] =  M[0] * omega[1] - M[1] * omega[0] + dr[0] * (-omega[0] * r[2] + omega[2] * r[0]) - dr[1] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * (-r[1] * gamma[0] + r[0] * gamma[1])

    H = dM[1] * gamma[0] + M[1] * dgamma[0] - dM[0] * gamma[1] - M[0] * dgamma[1]

    fx[0] = dM[0] / H
    fx[1] = dM[1] / H
    fx[2] = dM[2] / H
    fx[3] = dgamma[0] / H
    fx[4] = dgamma[1] / H
    fx[5] = dgamma[2] / H

def fix_integrals(alpha, beta, x):
    params[0] = alpha
    params[-2] = beta

    delta = params[0]
    J = np.array(params[1:4])
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    # Распакуем начальные условия
    gamma = x[2:]
    tmp_M = np.array([x[0]*(x[2]/x[3]),*x[:2]])

    # Посчитаем вектор r
    r = np.ndarray(3)
    r[0] = -params[4] * gamma[0] / gamma[2]
    r[1] = -params[5] * gamma[1] / gamma[2]
    r[2] = -params[6] + (params[4] * pow(gamma[0], 0.2e1) + params[5] * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.2e1) / 0.2e1

    # Посчитаем вектор omega
    JQ_kuz = np.ndarray((3,3))

    JQ_kuz[0][0] = pow(np.cos(delta), 0.2e1) * J[0] + pow(np.sin(delta), 0.2e1) * J[1] + pow(r[1], 0.2e1) + pow(r[2], 0.2e1)
    JQ_kuz[0][1] = -np.cos(delta) * J[0] * np.sin(delta) + np.sin(delta) * J[1] * np.cos(delta) - r[0] * r[1]
    JQ_kuz[0][2] = -r[0] * r[2]
    JQ_kuz[1][0] = -np.cos(delta) * J[0] * np.sin(delta) + np.sin(delta) * J[1] * np.cos(delta) - r[0] * r[1]
    JQ_kuz[1][1] = pow(np.sin(delta), 0.2e1) * J[0] + pow(np.cos(delta), 0.2e1) * J[1] + pow(r[0], 0.2e1) + pow(r[2], 0.2e1)
    JQ_kuz[1][2] = -r[1] * r[2]
    JQ_kuz[2][0] = -r[0] * r[2]
    JQ_kuz[2][1] = -r[1] * r[2]
    JQ_kuz[2][2] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) + J[2]

    JQ_rev = np.linalg.inv(JQ_kuz)

    omega = JQ_rev @ tmp_M

    # Отнормируем энергию
    M = np.zeros(3)
    for k in range(3):
        M[k] = tmp_M[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp_M, omega)))
    x[0] = M[1]
    x[1] = M[2]

    # Поправим геометрический интеграл
    tmp = gamma.copy()
    for k in range(3):
        x[k+2] = gamma[k]/(np.sqrt(np.dot(tmp, tmp)))

    # print(x)
    # En = (1/2) * np.dot(x[:3], omega) - g0 * np.dot(r, x[3:])
    # print("delta", params[0])
    # print("energy", En, "en_param", params[-2])
    # print("geom", np.dot(x[3:], x[3:]))
    return x

def test_F(x: np.ndarray, alpha: float, beta: float, step = 0.0025) -> np.ndarray:
    """
    Отображение Пуанкаре

    """
    dimension = 6

    params[0] = alpha
    params[-2] = beta

    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    ### Замена от переменных Андуайе-Депри
    l = x[0]
    LG = x[1]
    HG = x[2]

    m1 = (1-LG**2)**(1/2) * np.sin(l)
    m2 = (1-LG**2)**(1/2) * np.cos(l)
    m3 = LG
    M = np.array([m1,m2,m3]) # При такой замене вектор M будет единичной длины

    gamma1 = (HG * (1-LG**2)**(1/2) + LG*(1-HG**2)**(1/2)) * np.sin(l)
    gamma2 = (HG * (1-LG**2)**(1/2) + LG*(1-HG**2)**(1/2)) * np.cos(l)
    gamma3 = HG * LG - (1-LG**2)**(1/2) * (1-HG**2)**(1/2) * 1
    gamma = np.array([gamma1,gamma2,gamma3]) # гамма тоже единичный
    gamma = -gamma # надо повернуть

    QC = np.array([[np.cos(d), np.sin(d), 0],
                [-np.sin(d),np.cos(d),0],
                [0,0,1]])
    M = QC @ M
    gamma = QC @ gamma

    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    tmp = M.copy()
    for k in range(3):
        M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))
    fx = np.array([*M, *gamma])
    # print("После замены", fx)

    # Итерация отображения Пуанкаре
    # fx = np.array([x[0]*(x[2]/x[3]),*x])
    # print(fx)
    iter_num = 0
    while True:
        old_ps = fx[1] * fx[3] - fx[0] * fx[4]
        dverkStep(fx, dimension, diffFunc, params, step)
        new_ps = fx[1] * fx[3] - fx[0] * fx[4]
        iter_num += 1
        if iter_num > 1000:
            print("Poincare ERROR")
            return None
        if ((new_ps < 0) and (old_ps > 0) and iter_num!=1):
            dverkStep(fx, dimension, diffFuncP, params, -new_ps)
            break

    M = fx[:3]
    gamma = -fx[3:]

    QK = np.array([[np.cos(-d), np.sin(-d), 0],
                     [-np.sin(-d),np.cos(-d),0],
                     [0,0,1]])
    M = QK @ M
    gamma = QK @ gamma

    M = M / np.linalg.norm(M)

    L = M[2]
    G = np.sqrt(M[0]**2+M[1]**2+M[2]**2)
    H = np.dot(M,gamma)
    l = np.atan(M[0]/M[1]) # В этом случае пи надо прибавлять
    g = np.atan((M[1]*gamma[0]-M[0]*gamma[1])/(H*L/G-G*gamma[2]))

    # print("Результат итерации", l+np.pi, L/G, H/G)

    return np.array([l+np.pi, L/G, H/G])

def visualize_results(results: dict):
    """
    Визуализация результатов протягивания
    
    Parameters:
    -----------
    results : dict
        Результаты, возвращенные методом full_continuation
    """
    all_points = results['all_points']
    
    if not all_points:
        print("Нет точек для визуализации")
        return
    
    # Создаем массивы для графиков
    alphas = []
    betas = []
    iterations = []
    success_flags = []
    
    for point in all_points:
        alphas.append(point['alpha'])
        betas.append(point['beta'])
        iterations.append(point['iterations'])
        success_flags.append(1 if point['success'] else 0)
    
    alphas = np.array(alphas)
    betas = np.array(betas)
    iterations = np.array(iterations)
    success_flags = np.array(success_flags)
    
    # Создаем сетку для интерполяции
    if len(alphas) > 1 and len(betas) > 1:
        # Создаем регулярную сетку
        alpha_unique = np.unique(alphas)
        beta_unique = np.unique(betas)
        
        if len(alpha_unique) > 1 and len(beta_unique) > 1:
            alpha_grid, beta_grid = np.meshgrid(alpha_unique, beta_unique)
            iterations_grid = np.full_like(alpha_grid, np.nan, dtype=float)
            success_grid = np.full_like(alpha_grid, np.nan, dtype=float)
            
            # Заполняем сетку значениями
            point_dict = {}
            for point in all_points:
                key = (point['alpha'], point['beta'])
                point_dict[key] = point
            
            for i in range(alpha_grid.shape[0]):
                for j in range(alpha_grid.shape[1]):
                    key = (alpha_grid[i, j], beta_grid[i, j])
                    if key in point_dict:
                        iterations_grid[i, j] = point_dict[key]['iterations']
                        success_grid[i, j] = 1 if point_dict[key]['success'] else 0
    
    # Создаем графики
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # График 1: Распределение точек по параметрам
    ax = axes[0]
    # ax.set_xlim(-0.15,0.1)
    scatter1 = ax.scatter(alphas, betas, c=iterations, cmap='viridis', s=50, alpha=0.7)
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Beta')
    ax.set_title('Число итераций метода Ньютона')
    plt.colorbar(scatter1, ax=ax, label='Итерации')
    
    # График 2: Успешность сходимости
    ax = axes[1]
    # ax.set_xlim(-0.15,0.1)
    # Разделяем успешные и неуспешные точки
    success_mask = success_flags == 1
    if np.any(success_mask):
        ax.scatter(alphas[success_mask], betas[success_mask], 
                  c='green', s=50, alpha=0.7, label='Успех')
    if np.any(~success_mask):
        ax.scatter(alphas[~success_mask], betas[~success_mask], 
                  c='red', s=50, alpha=0.7, label='Неудача')
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Beta')
    ax.set_title('Успешность сходимости')
    ax.legend()
    
    # График 3: 3D график итераций
    if 'iterations_grid' in locals() and not np.all(np.isnan(iterations_grid)):
        ax = axes[2]
        # ax.set_xlim(-0.15,0.1)
        contour = ax.contourf(alpha_grid, beta_grid, iterations_grid, 
                             levels=20, cmap='viridis')
        ax.set_xlabel('Alpha')
        ax.set_ylabel('Beta')
        ax.set_title('Число итераций (интерполяция)')
        plt.colorbar(contour, ax=ax, label='Итерации')
    else:
        # Альтернативный график: гистограмма итераций
        ax = axes[2]
        ax.hist(iterations, bins=20, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Число итераций')
        ax.set_ylabel('Частота')
        ax.set_title('Распределение числа итераций')
    
    plt.tight_layout()
    plt.show()
    
    # Выводим статистику
    print("\nСтатистика:")
    print(f"Всего точек: {len(all_points)}")
    print(f"Успешных: {np.sum(success_flags)}")
    print(f"Неудачных: {len(all_points) - np.sum(success_flags)}")
    if len(iterations) > 0:
        print(f"Среднее число итераций: {np.mean(iterations):.2f}")
        print(f"Максимальное число итераций: {np.max(iterations)}")


# Пример использования
def main():
    # Параметры задачи
    dim = 3  # Фиксированная размерность

    # Начальные условия
    # x_start = np.array([-59.75329104, -38.86436549,  64.18350573, 0.27811237, 0.18088813, 0.9433626])
    # x_start = np.array([-38.86436549,  64.18350573, 0.27811237, 0.18088813, 0.9433626]) # M1 из уравнения секущей
    # x_start = np.array([3.59652, 0.667929, -0.3847012778])
    x_start = np.array([3.6507229058, 0.6691441589, -0.3847013043])

    # Параметры
    alpha0 = 0.485
    beta0 = 752
    
    # Диапазоны параметров
    # alpha_min, alpha_max = 0.41, 0.485
    # beta_min, beta_max = 735, 760
    alpha_min, alpha_max = 0.475, 0.485
    beta_min, beta_max = 750, 754
    
    # Шаги
    alpha_step = 0.0015
    beta_step = 0.5

    print(f"Результатом будет сетка {(alpha_max-alpha_min)/alpha_step}x{(beta_max-beta_min)/beta_step}")
    
    print("Инициализация системы...")
    print(f"Размерность: {dim}")
    print(f"Начальная точка: {x_start}")
    print(f"Начальные параметры: alpha={alpha0}, beta={beta0}")
    
    # Создаем экземпляр класса
    continuation = FixedPointContinuation(test_F, dim=dim)
    
    # Выполняем полное протягивание
    results = continuation.full_continuation(
        x_start=x_start,
        alpha0=alpha0, beta0=beta0,
        alpha_min=alpha_min, alpha_max=alpha_max, alpha_step=alpha_step,
        beta_min=beta_min, beta_max=beta_max, beta_step=beta_step,
        n_workers=8
    )
    # print(results.keys())
    # print(results['all_points_MG'])
    with open("points.txt", "w", encoding="utf-8") as f:
        for item in results['all_points_MG']:
            alpha = item['alpha']
            beta = item['beta']
            coords = item['x']

            # объединяем alpha, beta и координаты в одну строку
            line_values = [alpha, beta] + coords.tolist()

            # преобразуем в строку через запятую
            line = ", ".join(map(str, line_values))

            f.write(line + "\\n")
    # Визуализируем результаты
    visualize_results(results)
    
    # Пример доступа к результатам
    print("\nПримеры найденных точек:")
    for i, point in enumerate(results['all_points'][:3]):  # Показываем первые 3 точки
        print(f"Точка {i+1}: alpha={point['alpha']:.3f}, beta={point['beta']:.3f}, "
              f"итерации={point['iterations']}, успех={point['success']}")
        print(f"  x = {point['x']}")


if __name__ == "__main__":
    main()
    # x = np.array([-59.75329104, -38.86436549,  64.18350573, 0.27811237, 0.18088813, 0.9433626])
    # step = 0.0025
    # params = np.array([0.485, 2, 6, 7, 9, 4, 1, 752, 100])
    # delta_start = params[0]
    # E_start = params[-2]
    # res = []
    # for i in range(1000):
    #     # print("x1", x)
    #     res_point = test_F(x, delta_start, E_start, params, step)
    #     # print(res_point)
    #     res.append(res_point.copy())
    #     x = res_point.copy()
    #     # print("x", x)

    # res = np.array(res).T

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.scatter(res[2], res[3], res[4], s=1)
    # plt.show()
