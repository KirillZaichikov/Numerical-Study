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

fig1, ax1 = plt.subplots()

def main():
    mas_for_points = np.array([ [None] * dimension for i in range(int((integrate_time / step) / skip_for_phase))])
    mas_for_triger = np.array([ [None] * dimension for i in range(int(seq_len))])
    code = ""
    print("***********TIME for integrate:", (int(integrate_time / step)))

    start = time.time()
    i = 0
    calc_time = 0
    is_it_first = True
    while len(code) < seq_len:
        last_point = initial_point.copy()
        makeStep(initial_point, dimension, diffFunc, params, step)
        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = initial_point
        if last_point[2] < initial_point[2] and is_it_first:
            is_it_first = False
            mas_for_triger[len(code)] = initial_point
            if initial_point[0] > 0:
                code += "1"
            else:
                code += "0"
        elif last_point[2] >= initial_point[2]:
            is_it_first = True
        i += 1
        calc_time += step
    end = time.time()
    print(code)
    print("***********TIME:", end-start)

    mas_for_points = mas_for_points.T
    mas_for_points = mas_for_points[mas_for_points != None]
    mas_for_points = mas_for_points.reshape((dimension, int(len(mas_for_points) / dimension)))
 
    mas_for_triger = mas_for_triger.T
    mas_for_triger = mas_for_triger[mas_for_triger != None]
    mas_for_triger = mas_for_triger.reshape((dimension, int(len(mas_for_triger) / dimension)))

    # меняй индексы в зависимости от переменных
    if model_type == "map":
        ax.scatter(mas_for_points[0], mas_for_points[1], mas_for_points[2], s=1, c="black")
    else:
        ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linewidth=0.1)
    ax.scatter(mas_for_triger[0], mas_for_triger[1], mas_for_triger[2], s=5, c="red")
    ax1.plot([i for i in range(len(mas_for_triger[0]))], mas_for_triger[0], linewidth=0.3)

    print("Time for symbols: ", calc_time)
    plt.show()

if __name__ == "__main__":
    main()