import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
from params.params import *
from utils.integrator import *
from models.model import *


# for phase portret
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# менять тип графика в зависимости от размерности
# for poincare map
# fig3Dim = plt.figure()
# ax3Dim = fig3Dim.add_subplot(projection='3d')
fig3Dim, ax3Dim = plt.subplots()


def main():
    mas_for_poincare = np.array([ [None] * (dimension-1) for i in range(int(integrate_time / step))])
    mas_for_points = np.array([ [None] * dimension for i in range(int((integrate_time / step) / skip_for_phase))])
    
    print("***********TIME for integrate:", (int(integrate_time / step)))
    for2dCount = 0

    start = time.time()
    for i in range(int(skip_time / step)):
        makeStep(initial_point, dimension, diffFunc, params, step)
    print("***********END warm time:", time.time() - start)

    start = time.time()
    for i in range(int(integrate_time / step)):
        last_point = np.copy(initial_point)
        makeStep(initial_point, dimension, diffFunc, params, step)
        for2dCount = poincareMap(initial_point, last_point, mas_for_poincare, for2dCount)
        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = initial_point
    end = time.time()
    print("***********TIME:", end-start)

    start = time.time()
    mas_for_points, mas_for_poincare = mas_for_points.T, mas_for_poincare.T

    mas_for_poincare = mas_for_poincare[mas_for_poincare != None]
    mas_for_points = mas_for_points[mas_for_points != None]
    
    mas_for_poincare = mas_for_poincare.reshape((dimension-1, int(len(mas_for_poincare) / (dimension-1))))
    mas_for_points = mas_for_points.reshape((dimension, int(len(mas_for_points) / dimension)))
    print(mas_for_poincare.shape)
    print("***********FOR NUMPY:", time.time()-start)

    ax3Dim.scatter(mas_for_poincare[0], mas_for_poincare[1], s=1, c="black")
    # меняй индексы в зависимости от переменных
    if model_type == "map":
        ax.scatter(mas_for_points[0], mas_for_points[1], mas_for_points[2], s=1, c="black")
    else:
        ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linewidth=0.05)

    plt.show()


# менять в зависимости от размерности и переменной
def poincareMap(point, pre_point, massive, counter):
    if (point[2] < delta_params + eps) and (point[2] > delta_params - eps) and (pre_point[2] > delta_params+eps):
        massive[counter] = [point[0], point[1]] # или добавить еще переменную
        return counter + 1
    else:
        return counter


if __name__ == "__main__":
    main()