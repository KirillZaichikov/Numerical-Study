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
# fig, ax = plt.subplots()

def main(start_point):
    mas_for_points = np.array([ [None] * dimension for i in range(int((integrate_time / step) / skip_for_phase))])
    
    print("***********TIME for integrate:", (int(integrate_time / step)))

    start = time.time()
    for i in range(int(skip_time / step)):
        makeStep(start_point, dimension, diffFunc, params, step)
    print("***********END warm time:", time.time() - start)

    start = time.time()
    for i in range(int(integrate_time / step)):
        makeStep(start_point, dimension, diffFunc, params, step)
        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = start_point
    end = time.time()
    print("last point: ", start_point)
    print("***********TIME:", end-start)

    start = time.time()
    mas_for_points = mas_for_points.T
    mas_for_points = mas_for_points[mas_for_points != None]
    mas_for_points = mas_for_points.reshape((dimension, int(len(mas_for_points) / dimension)))
    print("***********FOR NUMPY:", time.time()-start)

    # меняй индексы в зависимости от переменных
    if model_type == "map":
        ax.scatter(mas_for_points[0], mas_for_points[1], mas_for_points[2], s=1, c="black")
    else:
        ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linewidth=0.2)


if __name__ == "__main__":
    start_point = np.array([initial_point])
    for i in start_point:
        main(i)
    plt.show()