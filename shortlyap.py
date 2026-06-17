import time
import math as m
from utils.linal import *
from params.params import *
from utils.integrator import *
import matplotlib.pyplot as plt
# from numba import jit

# @jit(nopython=True, cache=True)
def main_calc(initial_point, params):
    '''Вычислительный блок'''
    der = np.identity(dimension) * eps
    print(der)
    new_norm = np.zeros(lyap_num)

    # start = time.time()
    print("SKIP_TIME, STEP", skip_time, step)
    for i in range(int(skip_time/step)):
        makeStep(initial_point, dimension, diffFunc, params, step)
    # print("***********TIME for SKIP LYAPUNOV:", time.time() - start)

    # start = time.time()
    print("SKIP_VAR, STEP", skip_var_time, step)
    for i in range(int(skip_var_time/step)):
        for j in range(lyap_num):
            der[j] += initial_point
            makeStep(der[j], dimension, diffFunc, params, step)
            # print(der[0])
        makeStep(initial_point, dimension, diffFunc, params, step)
        # print(initial_point)
        for j in range(lyap_num):
            # if np.linalg.norm(der[j]) >= np.pi: 
            #     print("ERROR_WARM", der)
            der[j] -= initial_point
        ortVecs(der, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(der[j])
            der[j] = (der[j] / new_norm[j])*eps

    # batch_len = 10
    lyap_list = []
    for i in range(int(integrate_time/step/batch_len)):
        P = np.zeros(dimension, dtype=float)

        for iter in range(batch_len):
            for j in range(lyap_num):
                der[j] += initial_point
                makeStep(der[j], dimension, diffFunc, params, step)
            makeStep(initial_point, dimension, diffFunc, params, step)
            for j in range(lyap_num):
                # if new_norm[j] >= np.pi: print("ERROR_MAIN")
                der[j] -= initial_point
            ortVecs(der, dimension, lyap_num)
            for j in range(lyap_num):
                new_norm[j] = np.linalg.norm(der[j])
                der[j] = (der[j] / new_norm[j])*eps
                P[j] += m.log(new_norm[j]/eps)
        lyap_list.append((P/batch_len).copy())

    return lyap_list

def drawing_two_lines(y):
    '''КОД ДЛЯ ГРАФИКОВ'''
    # fig, ax = plt.subplots()
    # x = np.array([i for i in range(len(y[0]))])
    # clrs = ['red', 'magenta', 'blue']
    # for i, y1 in enumerate(y):
    #     print(len(x), len(y1))
    #     ax.plot(x, y1, color=clrs[i])
    # plt.show()

    '''КОД ДЛЯ ГИСТОГРАММ'''
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for line_num, lines_y in enumerate(y):
        for i in range(3):
            axes[line_num][i].hist(lines_y[i], bins=100, density=True)
            axes[line_num][i].set_title(rf'$\Lambda_{i+1}$', fontsize=16)
    
        axes[line_num][1].axvline(0, c='k', ls='--')
    plt.tight_layout()
    plt.savefig("hyst.png", dpi=400)
    plt.show()

def drawing_one_line(y):
    '''КОД ДЛЯ ГРАФИКОВ'''
    # fig, ax = plt.subplots()
    # x = np.array([i for i in range(len(y[0]))])
    # clrs = ['red', 'magenta', 'blue']
    # for i, y1 in enumerate(y):
    #     print(len(x), len(y1))
    #     ax.plot(x, y1, color=clrs[i])
    # plt.show()

    '''КОД ДЛЯ ГИСТОГРАММ'''
    fig, axes = plt.subplots(2, 3, figsize=(12, 4))
    for i in range(3):
        axes[i].hist(y[i], bins=100, density=True)
        axes[i].set_title(rf'$\Lambda_{i+1}$', fontsize=16)

    axes[1].axvline(0, c='k', ls='--')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    start_point = np.array(initial_point) # ЕСЛИ ХОЧЕШЬ В ДВЕ СТРОЧКИ, ТО И ТОЧЕК ТОЖЕ ДВЕ
    clrs = ["black", "red"]
    print(start_point)

    points = []
    for (i, j, k) in zip(start_point, clrs, params):
        print(i)
        lyap_l = main_calc(i, k) # СЮДА НУЖНО ПЕРЕДАВАТЬ ЕЩЕ И ПАРАМЕТРЫ
        points.append(np.array(lyap_l).T)
    print(points[0])
    # print(points[1])
    drawing_two_lines(points)

    # ОБЫЧНЫЕ В ОДНУ СТРОЧКУ
    # for (i, j) in zip(start_point, clrs):
    #     print(i)
    #     lyap_l = main_calc(i)
    #     arr = np.array(lyap_l)
    # drawing(arr)