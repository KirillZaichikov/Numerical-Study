import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
from params.params import *
from utils.integrator import *
from models.model import *


# for phase portret
# fig = plt.figure(figsize=(8,8))
# ax = fig.add_subplot(projection='3d')
# scale = 0.02
# ax.set_xlim(-scale*0.1, scale)
# ax.set_ylim(-scale, scale*0.1)
# ax.set_zlim(-np.pi, np.pi)


# fig, ax = plt.subplots()

fig = plt.figure(figsize=(6, 6))
ax2d_2 = fig.add_subplot()
ax2d_2.tick_params(axis='both', labelsize=15)

# Создаем сетку 2×2
# gs = fig.add_gridspec(2, 2)

# # 3D-график (занимает всю левую колонку)
# ax3d = fig.add_subplot(gs[:, 0], projection='3d')
# ax3d.get_xaxis().set_visible(False)
# ax3d.get_yaxis().set_visible(False)
# ax3d.get_zaxis().set_visible(False)
# ax3d.axis('off')
# ax3d.grid(True)

# 2D-графики справа
# ax2d_1 = fig.add_subplot(gs[0, 1])
# ax2d_2 = fig.add_subplot(gs[1, 1])

# ticks = [
#     -2*np.pi, -3*np.pi/2, -np.pi, -np.pi/2, 
#     0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi
# ]
# labels = [
#     r'$-2\pi$',
#     r'$-\frac{3\pi}{2}$',
#     r'$-\pi$',
#     r'$-\frac{\pi}{2}$',
#     r'$0$',
#     r'$\frac{\pi}{2}$',
#     r'$\pi$',
#     r'$\frac{3\pi}{2}$',
#     r'$2\pi$'
# ]
# ax2d_2.set_xticks(ticks)
# ax2d_2.set_xticklabels(labels, fontsize=16)

# ax2d_2.set_xlim(-np.pi, np.pi)
# ax2d_2.set_xlabel("teta")
# ax2d_2.set_ylabel("ksi")

def drawing(mas_for_points, clr):
    # меняй индексы в зависимости от переменных
    if model_type == "map":
        # ax.plot(mas_for_points[0][::10], mas_for_points[1][::10], mas_for_points[2][::10], linestyle="", marker="o", markersize=0.3, color="black")
        # ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linestyle="", marker="o", markersize=0.3, color="black", rasterized=True)
        # ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linestyle="", marker="o", markersize=2, color="black", rasterized=True)

        ax2d_2.plot(mas_for_points[2], mas_for_points[0], linestyle="", marker="o", markersize=0.1, color=clr, rasterized=True)
        
        # ax2d_2.plot(mas_for_points[2], mas_for_points[0], linestyle="", marker="o", markersize=2, color="black", rasterized=True)
        # ax2d_2.plot(mas_for_points[2], mas_for_points[0], ',k', alpha=0.25) #0.25
        # ax2d_2.scatter(mas_for_points[2], mas_for_points[0], marker=',', color=clr, alpha=0.25) #0.25
        # ax2d_2.plot(mas_for_points[2], mas_for_points[0], linestyle="", marker=',', color=clr, alpha=0.25)
        # plt.tight_layout()
        plt.savefig('high_res_plot.png', dpi=400, bbox_inches="tight")

        # ax.plot(mas_for_points[0], mas_for_points[1], linestyle="", marker="o", markersize=1, color="black")
    else:
        ax.plot(mas_for_points[0], mas_for_points[1], mas_for_points[2], linewidth=0.2, rasterized=True)


def main(start_point, clr):
    mas_for_points = np.array([ [None] * dimension for i in range(int((integrate_time / step) / skip_for_phase))])
    
    print("***********TIME for integrate:", (int(integrate_time / step)))

    start = time.time()
    for i in range(int(skip_time / step)):
        makeStep(start_point, dimension, diffFunc, params, step)
    print("***********END warm time:", time.time() - start)

    start = time.time()
    # print(start_point)
    for i in range(int(integrate_time / step)):
        makeStep(start_point, dimension, diffFunc, params, step)
        # print(start_point)
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

    drawing(mas_for_points, clr)


if __name__ == "__main__":
    start_point = np.array(initial_point)
    clrs = ["black", "red"]
    print(start_point)
    for (i, j) in zip(start_point, clrs):
        print(i)
        main(i, j)
    print("*********** calc is ending")
    plt.savefig("plot.pdf", dpi=300)
    plt.show()