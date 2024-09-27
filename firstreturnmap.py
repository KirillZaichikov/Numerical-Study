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

# Менять в зависимости от размерности системы
# for poincare map
# fig3Dim = plt.figure()
# ax3Dim = fig3Dim.add_subplot(projection='3d')
figForMap, axForMap = plt.subplots()

INDEX = 0 # какую по счету переменную запоминаем


def main():
    # 3.6 for shimizu
    alpha = 2 * np.pi #/ 5.8

    rotate_matrix = np.array([[np.cos(alpha), np.sin(alpha)],
                              [-np.sin(alpha), np.cos(alpha)]])
    print("***********ROTATE MATRIX")
    print(rotate_matrix)

    # mas_for1D = np.ndarray((2000, 2))
    mas_for1D = np.array([ [None] * 2 for i in range(2000)])
    mas_for_points = np.array([ [None] * dimension for i in range(int((integrate_time / step) / skip_for_phase))])
    
    print("***********TIME for integrate:", (int(integrate_time / step)))
    for1dCount = 0
    domestic_var = None
    max_step = 0

    start = time.time()
    for i in range(int(skip_time / step)):
        makeStep(initial_point, dimension, diffFunc, params, step)
    for i in range(int(integrate_time / step)):
        last_point = np.copy(initial_point)
        makeStep(initial_point, dimension, diffFunc, params, step)
        if abs(last_point[2] - initial_point[2]) > max_step:
            max_step = abs(last_point[2] - initial_point[2])
        for1dCount, domestic_var = poincareMap(initial_point, last_point, for1dCount, mas_for1D, domestic_var)
        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = initial_point
    print(max_step)
    print("***********TIME:", time.time()-start)

    start = time.time()

    mas_for1D = mas_for1D.T
    mas_for1D = mas_for1D[mas_for1D != None]
    mas_for1D = mas_for1D.reshape((2, int(len(mas_for1D) / 2)))

    mas_for_points = mas_for_points.T
    mas_for_points = mas_for_points[mas_for_points != None]
    mas_for_points = mas_for_points.reshape((dimension, int(len(mas_for_points) / dimension)))
    print("***********FOR NUMPY:", time.time()-start)
    # mas_for_property = mas_for_property[:for2dCount]
    # mas_for_propertyT = mas_for_property.T
    # mas_for_propertyT = rotate_matrix @ mas_for_propertyT
    # mas_for_poincare = mas_for_poincare[mas_for_poincare != 0
    axForMap.scatter(mas_for1D[0], mas_for1D[1], s=1)
    if model_type == "map":
        ax.scatter(mas_for_points[0], mas_for_points[1], mas_for_points[2], s=1, c="black")
    else:
        ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linewidth=0.1)

    plt.show()


def poincareMap(point, pre_point, counter, mas_for1D, domestic_var):
    if (point[2] < delta_params + eps) and (point[2] > delta_params - eps) and (pre_point[2] < delta_params - eps):
        if domestic_var == None:
            domestic_var = point[INDEX]
            return counter, domestic_var
        else:
            if domestic_var > 0:
                if point[INDEX] > 0:
                    mas_for1D[counter] = [domestic_var, point[INDEX]]
                else:
                    mas_for1D[counter] = [domestic_var, -point[INDEX]]
                domestic_var = point[INDEX]
            else:
                domestic_var = point[INDEX]
                return counter, domestic_var

        return counter + 1, domestic_var
    else:
        return counter, domestic_var


if __name__ == "__main__":
    main()