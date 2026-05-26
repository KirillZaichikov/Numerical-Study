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
    while start_point[0] > 2*np.pi:
        start_point[0] -= 2*np.pi
    while start_point[0] <0:
        start_point[0] += 2*np.pi

    while start_point[1] > 2*np.pi:
        start_point[1] -= 2*np.pi
    while start_point[1] < 0:
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
def get_first_cross(start_point,mas_for_points):

    flag = 0
    for i in range(int(integrate_time / step)):
        old_point = start_point.copy()
        makeStep(start_point, dimension_after_replace, diffFunc, params, step)
        norm_solution(start_point)

        # if (old_point[1] > crossection and start_point[1] < crossection) or (old_point[1] < crossection and start_point[1] > crossection):
        # if (old_point[1] < crossection and start_point[1] > crossection):
        if (old_point[1] > crossection and start_point[1] < crossection and abs(start_point[1]-old_point[1])<np.pi) or \
            (old_point[1] < crossection and start_point[1] > crossection and abs(start_point[1]-old_point[1])<np.pi):
            flag = 1
            H = - get_omega(start_point, params) * ctg(start_point[2]) * np.sin(start_point[0])
            # print("FIND CROSS")
            # print("H", -(start_point[1] - np.pi))
            # print("old", old_point)
            # print("new", start_point)
            makeStep(start_point, dimension_after_replace, diffPoincare, params, -(start_point[1] - crossection), H)
            # print(i, start_point)
            break
        
        if np.isnan(start_point[0]):
            # print("old_point",old_point)
            break

        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = start_point

    # print("last point: ", start_point)
    if flag:
        return start_point.copy()
    else:
        return np.array([np.nan,np.nan,np.nan])
        

figPoin = plt.figure(figsize=(6, 6))
axPoin = figPoin.add_subplot()
axPoin.set_xlim(0,2*np.pi)
axPoin.set_ylim(0,2*np.pi)

if __name__ == "__main__":
    start = time.time()
    phi = np.linspace(eps, 2*np.pi-eps, DISCR)
    arr = np.empty((len(phi), int((integrate_time / step) / skip_for_phase), dimension_after_replace))
    
    result_list = []
    for number, i in enumerate(phi):
        result_list.append(np.array([eps, i, eps]))

    start_point = np.array(result_list)
    mas_for_points = arr
    
    cross_points = []
    for number, i in enumerate(start_point):
        print(number)
        cross_points.append(get_first_cross(i, mas_for_points[number]))

    # for number, elem in enumerate(mas_for_points):
    #     print(number)
    #     nelem = elem.T
    #     nelem = nelem[nelem != None]
    #     nelem = nelem.reshape((dimension, int(len(nelem) / dimension)))
    #     ax.plot(nelem[0], nelem[1], nelem[2], linestyle="", marker="o", markersize=0.1, color="black")

    cross_points = np.array(cross_points)
    print("колво точек на сечении", len(cross_points))
    cross_points = cross_points[~np.isnan(cross_points).any(axis=1)].T
    # print(cross_points)
    axPoin.plot(cross_points[0], cross_points[2], linestyle="", marker="o", markersize=0.5, color="black")

    plotter = pv.Plotter()

    plotter.show_bounds(bounds=(0, 2*np.pi,   # X
                                0, 2*np.pi,   # Y
                                0, np.pi),     # Z
                        grid='back',        # сетка
                        location='outer',   # снаружи
                        all_edges=True      # все рёбра
    )

    colors = plt.cm.plasma(np.linspace(0, 1, len(mas_for_points)))
    points_number = 0
    for i, (traj, color) in enumerate(zip(mas_for_points, colors)):
        # if i == 0:
            points_number += len(traj)
            points = pv.PolyData(traj)
            plotter.add_mesh(points, color=color[:3], point_size=1)
    print(points_number)
    print('result time', time.time()-start)
    plotter.show()

    plt.show()