import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
import pyvista as pv
from numba import jit

from params import *
from integrator import *
from system import *

np.set_printoptions(suppress=True, precision=10)

@jit(nopython=True, cache=True)
def norm_solution(start_point):
    if start_point[0] > 2*np.pi:
        start_point[0] = start_point[0] % (np.pi * 2)
    elif start_point[0] <0:
        start_point[0] += 2*np.pi

    if start_point[1] > 2*np.pi:
        start_point[1] = start_point[1] % (np.pi * 2)
    elif start_point[1] < 0:
        start_point[1] += 2*np.pi

@jit(nopython=True, cache=True)
def get_omega(start_point, params):
    ro, mpar, i1, i2, a0, g, eps = params
    ksi, phi, teta = start_point
    omega2 = 2*(eps - mpar * g * m.sin(teta) * (ro + a0*m.sin(phi)) ) / \
        (i1 * m.cos(ksi-phi) ** 2 + i2 * m.sin(ksi-phi) ** 2 + \
            mpar * m.cos(teta)**2*(ro*m.cos(ksi)-a0 * m.sin(ksi-phi))**2)
    # print(omega2)
    omega = m.sqrt(omega2)
    return omega

@jit(nopython=True, cache=True)
def get_cross(start_point):

    flag = 0
    for i in range(int(integrate_time / step)):
        old_point = start_point.copy()
        makeStep(start_point, dimension, diffFunc, params, step)
        norm_solution(start_point)

        # if (old_point[1] > crossection and start_point[1] < crossection) or (old_point[1] < crossection and start_point[1] > crossection):
        # if (old_point[1] < crossection and start_point[1] > crossection):
        if (old_point[1] > crossection and start_point[1] < crossection and abs(start_point[1]-old_point[1])<np.pi) or \
            (old_point[1] < crossection and start_point[1] > crossection and abs(start_point[1]-old_point[1])<np.pi):
            flag = 1
            H = - get_omega(start_point, params) * ctg(start_point[2]) * m.sin(start_point[0])
            makeStep(start_point, dimension, diffPoincare, params, -(start_point[1] - np.pi), H)
            break
        
        if np.isnan(start_point[0]):
            # print("old_point",old_point)
            break

    # print("last point: ", start_point)
    if flag:
        return start_point.copy()
    else:
        return np.array([np.nan,np.nan,np.nan])
        

figPoin = plt.figure(figsize=(6, 6))
axPoin = figPoin.add_subplot()
axPoin.set_xlim(0, 2*np.pi)
axPoin.set_ylim(0, 1)
axPoin.set_xlabel("dzeta")
axPoin.set_ylabel("teta")

if __name__ == "__main__":
    start = time.time()
    phi = np.linspace(eps, 2*np.pi-eps, DISCR)
    arr = np.empty((len(phi), int((integrate_time / step) / skip_for_phase), dimension))
    
    result_list = []
    for number, i in enumerate(phi):
        result_list.append(np.array([eps, i, eps]))

    start_point = np.array(result_list)
    mas_for_points = arr
    
    cross_points = []
    for number, i in enumerate(start_point):
        print(number)
        start = i
        for j in range(iter_num):
            start = get_cross(start)
            cross_points.append(start.copy())

    cross_points = np.array(cross_points)
    print("колво точек на сечении", len(cross_points))
    cross_points = cross_points[~np.isnan(cross_points).any(axis=1)].T
    # print(cross_points)
    axPoin.plot(cross_points[0], cross_points[2], linestyle="", marker="o", markersize=0.5, color="black")

    # plotter = pv.Plotter()

    # plotter.show_bounds(bounds=(0, 2*np.pi,   # X
    #                             0, 2*np.pi,   # Y
    #                             0, np.pi),     # Z
    #                     grid='back',        # сетка
    #                     location='outer',   # снаружи
    #                     all_edges=True      # все рёбра
    # )

    # colors = plt.cm.plasma(np.linspace(0, 1, len(mas_for_points)))
    # points_number = 0
    # for i, (traj, color) in enumerate(zip(mas_for_points, colors)):
    #     # if i == 0:
    #         points_number += len(traj)
    #         points = pv.PolyData(traj)
    #         plotter.add_mesh(points, color=color[:3], point_size=1)
    # print(points_number)
    # print('result time', time.time()-start)
    # plotter.show()

    plt.show()