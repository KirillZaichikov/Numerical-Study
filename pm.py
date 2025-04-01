import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
from params.params import *
from utils.integrator import *
from models.model import *


step = 0.01
params = np.array([1.3131249999999999, 0.35854032258064517, 0.2], dtype=np.longdouble)
integrate_time = 10000
skip_for_phase = 1

def main():
    # coord_list_for_poincare = np.array([ [None] * (dimension-1) for i in range(int(integrate_time / step / skip_for_phase))], dtype=np.longdouble)
    # coord_list_for_phase = np.array([ [None] * dimension for i in range(int((integrate_time / step) / skip_for_phase))], dtype=np.longdouble)
    coord_list_for_poincare = []
    coord_list_for_phase = []

    start = time.time()
    for i in range(int(skip_time / step)):
        makeStep(initial_point, dimension, diffFunc, params, step)
    print("***********END warm time:", time.time() - start)

    start = time.time()
    for i in range(int(integrate_time/step)):
        z_last = initial_point[2] - crossection
        dverkStep(initial_point, dimension, ShimizuX3_3D_flow, params, step)
        if (initial_point[2] - crossection > 0) and (z_last < 0):
            H = - params[1] * initial_point[2] + initial_point[0] ** 2
            dverkStep(initial_point, dimension, ShimizuX3_3D_flow, params, -(initial_point[2] - crossection), H)
            coord_list_for_poincare.append(list(initial_point[:2]))
        # coord_list_for_phase.append(list(initial_point))
    print("***********TIME:", time.time() - start)

    # подготовка массивов
    transpose_list_for_phase = (np.array(coord_list_for_phase, dtype=np.longdouble)).T
    transpose_points_for_poincare = (np.array(coord_list_for_poincare, dtype=np.longdouble)).T
    np.save('2dimmaps/'+'a_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:9], transpose_points_for_poincare)
    
    # отрисовка фазового
    # fig_phase = plt.figure()
    # ax_phase = fig_phase.add_subplot(projection='3d')
    # ax_phase.plot(transpose_list_for_phase[0][::10], transpose_list_for_phase[1][::10], transpose_list_for_phase[2][::10]) # отрисовка фазового портрета
    # ax_phase.scatter(transpose_points_for_poincare[0], transpose_points_for_poincare[1], crossection, s=5, c="red") # отрисовка точек сечения которые идут для отображения первого возвращения

    # отрисовка сечения
    fig_poincare, ax_poincare = plt.subplots()
    ax_poincare.scatter(transpose_points_for_poincare[0], transpose_points_for_poincare[1], s=5, c="black") # отрисовка точек сечения
    plt.show()

# def Lorenz_4D_flow(state, res, params):  # Lorenz4D
#     H = (- params[2] * state[2] + params[3] * state[3] + state[0] * state[1])
#     res[0] = (params[0] * ( - state[0] + state[1] ) ) / H
#     res[1] = (state[0] * ( params[1] - state[2] ) - state[1]) / H
#     res[2] = (- params[2] * state[2] + params[3] * state[3] + state[0] * state[1]) / H
#     res[3] = (- params[2] * state[3] - params[3] * state[2]) / H

def ShimizuX3_3D_flow(state, res, params, H) -> None:
    # Param - alpha lamda B 
    res[0] = state[1] / H
    res[1] = (params[2] * pow ( state[0] , 3.0 ) - params[0] * state[1] - state[0] * state[2] + state[0]) / H
    res[2] = (- params[1] * state[2] + pow ( state[0] , 2.0 )) / H

def dverkStep(val, dimension, diffFunc, params, step, H=1) -> None:
    k1, k2, k3, k4, k5, k6, k7, k8 = np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), \
        np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension)
    arg = np.zeros(dimension)
    diffFunc(val, k1, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step / 6. * k1[j]
    
    diffFunc(arg, k2, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j])
    
    diffFunc(arg, k3, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j])
    
    diffFunc(arg, k4, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j])
    
    diffFunc(arg, k5, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j])
    
    diffFunc(arg, k6, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j])
    
    diffFunc(arg, k7, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j])
    
    diffFunc(arg, k8, params, H)
    for j in range(dimension):
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j])

if __name__ == "__main__":
    main()