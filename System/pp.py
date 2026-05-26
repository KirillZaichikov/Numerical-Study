import numpy as np
import math as m
import matplotlib.pyplot as plt
import time

from params import *
from integrator import *
from system import *


np.set_printoptions(suppress=True, precision=10)

# fig = plt.figure(figsize=(6, 6))
# ax = fig.add_subplot()
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim(0, 2*np.pi)
ax.set_ylim(0, 2*np.pi)
ax.set_zlim(0, np.pi)

ax.set_xlabel("ksi")
ax.set_ylabel("phi")
ax.set_zlabel("teta")

def norm_solution(start_point):
    while start_point[0] > 2*np.pi:
        start_point[0] -= 2*np.pi
    while start_point[0] < 0:
        start_point[0] += 2*np.pi
    while start_point[1] > 2*np.pi:
        start_point[1] -= 2*np.pi
    while start_point[1] < 0:
        start_point[1] += 2*np.pi
    # while start_point[2] > np.pi:
    #     start_point[2] -= np.pi
    # while start_point[2] < 0:
    #     start_point[2] += np.pi

def main(start_point):
    mas_for_points = np.array([ [None] * dimension for i in range(int((integrate_time / step) / skip_for_phase))])
    
    print("***********TIME for integrate:", (int(integrate_time / step)))

    start = time.time()
    for i in range(int(skip_time / step)):
        makeStep(start_point, dimension, diffFunc, params, step)
        # start_point[:2] = start_point[:2] % (2 * np.pi)
        norm_solution(start_point)
        if np.isnan(start_point[0]):
            break
        # print(start_point)
    print("***********END warm time:", time.time() - start)

    start = time.time()
    for i in range(int(integrate_time / step)):
        makeStep(start_point, dimension, diffFunc, params, step)
        norm_solution(start_point)
        # start_point[:2] = start_point[:2] % (2 * np.pi)
        # print(start_point)
        if np.isnan(start_point[0]):
            break
        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = start_point
    end = time.time()
    print(mas_for_points[-10:])
    print("last point: ", start_point)
    print("***********TIME:", end-start)

    start = time.time()
    mas_for_points = mas_for_points.T
    mas_for_points = mas_for_points[mas_for_points != None]
    mas_for_points = mas_for_points.reshape((dimension, int(len(mas_for_points) / dimension)))

    # print(mas_for_points[2])
    # меняй индексы в зависимости от переменных
    if model_type == "map":
        # ax3d.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linestyle="", marker="o", markersize=0.3, color="black", rasterized=True)

        ax2d_2.plot(mas_for_points[2], mas_for_points[0], linestyle="", marker="o", markersize=0.05, color="black", rasterized=True)
        
        # ax2d_2.plot(mas_for_points[2], mas_for_points[0], linestyle="", marker="o", markersize=2, color="black", rasterized=True)
        # ax2d_2.plot(mas_for_points[2], mas_for_points[0], ',k', alpha=0.25) #0.25
        plt.savefig('high_res_plot.png', dpi=400)
        # ax.plot(mas_for_points[0], mas_for_points[1], linestyle="", marker="o", markersize=1, color="black")
    else:
        # ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linewidth=0.2, rasterized=True)
        ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linestyle="", marker="o", markersize=0.1, color="black")


if __name__ == "__main__":
    start_point = np.array([initial_point])
    print("start_point", start_point)
    for i in start_point:
        main(i)
    print("*********** calc is ending")
    plt.savefig("plot.pdf", dpi=300)
    plt.show()