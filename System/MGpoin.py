import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
import pyvista as pv

from params import *
from integrator import *
from system import *


np.set_printoptions(suppress=True, precision=10)

def norm_solution(start_point):
    if start_point[0] > 2*np.pi:
        start_point[0] -= 2*np.pi
    elif start_point[0] <0:
        start_point[0] += 2*np.pi

    if start_point[1] > 2*np.pi:
        start_point[1] -= 2*np.pi
    elif start_point[1] < 0:
        start_point[1] += 2*np.pi

def replace_var(start_point):
    '''
    MG -> dzeta, phi, teta
    '''
    omega1, omega2, _g1, _g2, _g3 = start_point
    _phi = np.atan2(_g1, _g2)
    _teta = np.acos(_g3)
    _ksi = np.atan2(omega2, omega1)
    point = np.array([_ksi+_phi,_phi,_teta])
    norm_solution(point)
    return point.copy()

def calc_energy_dpt(start_point,params):
    dzeta, phi, teta = start_point
    ksi = dzeta-phi
    ro, mpar, i1, i2, a0, g, en = params

    omega = np.sqrt(2*(en - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
        (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
            mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2))
    
    E = ((omega ** 2) / 2) * (i1 * np.cos(ksi)**2 + i2 * np.sin(ksi)**2) + \
         (mpar / 2) * np.cos(teta)**2 * omega**2 * (ro * np.cos(phi+ksi) - \
                                                    -a0*np.sin(ksi)) ** 2 + \
         mpar * g * np.sin(teta) * (ro + a0 * np.sin(phi))
    # print("E(dpt)", E)

def back_replace_var(start_point, params):
    '''
    dzeta,phi,teta -> MG
    '''
    dzeta, phi, teta = start_point
    ro, mpar, i1, i2, a0, g, en = params

    gamma1=np.sin(teta)*np.sin(phi)
    gamma2=np.sin(teta)*np.cos(phi)
    gamma3=np.cos(teta)

    omega = np.sqrt(2*(en - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
        (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
            mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2))

    omega1 = omega * np.cos(dzeta-phi)
    omega2 = omega * np.sin(dzeta-phi)

    return np.array([omega1, omega2, gamma1, gamma2, gamma3]).copy()

def calc_energy(start_point, params):
    ro, mpar, i1, i2, a0, g, en = params
    omega = np.array([start_point[0], start_point[1], 0])
    gamma = start_point[2:]
    r = np.array([-(ro*gamma[0])/(np.sqrt(1-gamma[2]**2))-a0,
                  -(ro*gamma[1])/(np.sqrt(1-gamma[2]**2)),
                  0])

    J = np.eye(3) + mpar * np.outer((np.cross(r, gamma)), (np.cross(r, gamma)))
    E = (1/2) * np.dot(omega, J @ omega) - mpar*g*np.dot(r,gamma)
    # print("Energy", E)
    return E

def calc_geom_integral(start_point):
    geom = np.linalg.norm(start_point[2:])
    # print("geom", geom)
    return geom

def get_omega(start_point, params):
    ro, mpar, i1, i2, a0, g, eps = params
    ksi, phi, teta = start_point
    omega2 = 2*(eps - mpar * g * m.sin(teta) * (ro + a0*m.sin(phi)) ) / \
        (i1 * m.cos(ksi-phi) ** 2 + i2 * m.sin(ksi-phi) ** 2 + \
            mpar * m.cos(teta)**2*(ro*m.cos(ksi)-a0 * m.sin(ksi-phi))**2)
    # print(omega2)
    omega = m.sqrt(omega2)
    return omega

def get_first_cross(start_point, mas_for_points):
    start = time.time()
    dpt = start_point
    for i in range(int(integrate_time / step)):
        old_point = start_point.copy()
        start_point = back_replace_var(dpt, params)
        makeStep(start_point, dimension, diffFunc, params, step)
        dpt = replace_var(start_point)
        norm_solution(start_point)

        if (old_point[1] > crossection and dpt[1] < crossection and abs(dpt[1]-old_point[1])<np.pi) or \
           (old_point[1] < crossection and dpt[1] > crossection and abs(dpt[1]-old_point[1])<np.pi):
            H = - get_omega(dpt, params) * ctg(dpt[2]) * np.sin(dpt[0])
            print("FIND CROSS")
            makeStep(dpt, dimension_after_replace, diffPoincare, params, -(dpt[1] - np.pi), H)
            print(i, dpt)
            break

        if np.isnan(start_point[0]):
            print("NANANANAN")
            break
        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = dpt

    end = time.time()
    print("***********TIME:", end-start)
    return dpt.copy()

figPoin = plt.figure(figsize=(6, 6))
axPoin = figPoin.add_subplot()
axPoin.set_xlim(0,2*np.pi)
axPoin.set_ylim(0,2*np.pi)

if __name__ == "__main__":

    phi = np.linspace(eps, 2*np.pi-eps, DISCR)
    mas_for_points = np.empty((len(phi), int((integrate_time / step) / skip_for_phase), dimension_after_replace))
    result_list = []
    for number, i in enumerate(phi):
        result_list.append(np.array([eps, i, eps]))

    start_points = np.array(result_list)
    
    cross_points = []
    for number, i in enumerate(start_points):
        print(number)
        cross_points.append(get_first_cross(i, mas_for_points[number]))
    print(cross_points)
    # for number, elem in enumerate(mas_for_points):
    #     print(number)
    #     nelem = elem.T
    #     nelem = nelem[nelem != None]
    #     nelem = nelem.reshape((dimension, int(len(nelem) / dimension)))
    #     ax.plot(nelem[0], nelem[1], nelem[2], linestyle="", marker="o", markersize=0.1, color="black")

    cross_points = np.array(cross_points).T
    # cross_points = cross_points[cross_points!=None]
    # cross_points = cross_points[cross_points!=np.nan]
    # cross_points = cross_points.reshape((dimension_after_replace, int(len(cross_points) / dimension_after_replace)))
    print(cross_points)
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
            points_number += len(traj)
            points = pv.PolyData(traj)
            plotter.add_mesh(points, color=color[:3], point_size=1)
    print(points_number)
    plotter.show()

    plt.show()