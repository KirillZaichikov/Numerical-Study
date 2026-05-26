import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
import pyvista as pv

from params import *
from integrator import *
from system import *


np.set_printoptions(suppress=True, precision=10)

# fig = plt.figure(figsize=(6, 6))
# ax = fig.add_subplot()
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim(0, 2*np.pi)
# ax.set_ylim(0, 2*np.pi)
# ax.set_zlim(0, np.pi)

ax.set_xlabel("phi")
ax.set_ylabel("omega1")
ax.set_zlabel("omega2")

@jit(nopython=True, cache=True)
def norm_solution(start_point):
    if start_point[2] > 2*np.pi:
        start_point[2] -= 2*np.pi
    elif start_point[2] < 0:
        start_point[2] += 2*np.pi

@jit(nopython=True, cache=True)
def get_omega(start_point, params):
    ro, mpar, i1, i2, a0, g, eps = params
    dzeta, phi, teta = start_point
    omega2 = 2*(eps - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
        (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
            mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2)
    # print(omega2)
    omega = m.sqrt(omega2)
    return omega

@jit(nopython=True, cache=True)
def replace_var(start_point):
    '''
    MG -> omega1, omega2, phi
    '''
    omega1, omega2, _g1, _g2, _g3 = start_point
    # _ksi = np.atan2(omega2, omega1)
    _phi = np.atan2(_g1, _g2)
    # dzeta = _ksi - _phi
    point = np.array([omega1, omega2, _phi])
    norm_solution(point)
    return point.copy()

@jit(nopython=True, cache=True)
def ood_to_oop(start_point):
    '''
    omega1, omega2, dzeta -> omega1, omega2, phi
    '''
    omega1, omega2, dzeta = start_point
    _ksi = np.atan2(omega2, omega1)
    # _phi = np.atan2(_g1, _g2)
    _phi = _ksi - dzeta
    # dzeta = _ksi - _phi
    point = np.array([omega1, omega2, _phi])
    norm_solution(point)
    return point.copy()

@jit(nopython=True, cache=True)
def oop_to_ood(start_point):
    '''
    omega1, omega2, phi -> omega1, omega2, dzeta
    '''
    omega1, omega2, _phi = start_point
    _ksi = np.atan2(omega2, omega1)
    # _phi = np.atan2(_g1, _g2)
    # _phi = _ksi - dzeta
    dzeta = _ksi - _phi
    point = np.array([omega1, omega2, dzeta])
    # norm_solution(point)
    return point.copy()

@jit(nopython=True, cache=True)
def back_replace_var(start_point, params):
    '''
    omega1,omega2,phi -> OmegaGamma
    '''
    omega1, omega2, phi = start_point
    ro, mpar, i1, i2, a0, g, en = params

    ksi = np.atan2(omega2, omega1)
    # phi = dzeta - ksi
    dzeta = ksi + phi
    omega = np.sqrt(omega1 ** 2 + omega2 ** 2)

    # Результат зависит от начального условия
    # teta0 = np.pi-eps
    teta0 = 0
    for i in range(10):
        fx0 = ((omega ** 2) / 2) * (i1*np.cos(ksi)**2+i2*np.sin(ksi)**2) + \
            (mpar/2) * np.cos(teta0) ** 2 *(omega**2) * (ro * np.cos(phi+ksi) - \
                                                        a0*np.sin(ksi))**2 +\
            mpar*g*np.sin(teta0) * (ro + a0*np.sin(phi)) - en
        dteta = teta0 + 1e-7
        fdx0 = ((((omega ** 2) / 2) * (i1*np.cos(ksi)**2+i2*np.sin(ksi)**2) + \
            (mpar/2) * np.cos(dteta) ** 2 *(omega**2) * (ro * np.cos(phi+ksi) - \
                                                        a0*np.sin(ksi))**2 +\
            mpar*g*np.sin(dteta) * (ro + a0*np.sin(phi)) - en) - fx0) / 1e-7
        teta1 = teta0 - fx0/fdx0
        teta0 = teta1

    # print(fx0)
    gamma1=np.sin(teta0)*np.sin(phi)
    gamma2=np.sin(teta0)*np.cos(phi)
    gamma3=np.cos(teta0)

    return np.array([omega1, omega2, gamma1, gamma2, gamma3]).copy()

@jit(nopython=True, cache=True)
def replace_dpt(start_point, params):
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

@jit(nopython=True, cache=True)
def get_first_cross(start_point, mas_for_points, count_m, n):
    '''
    start_point: [omega1, omega2, dzeta]
    '''
    # start = time.time()

    flag = 0
    ood = start_point.copy()

    oop = ood_to_oop(ood)
    start_point = back_replace_var(oop, params)
    # dpt = replace_dpt(start_point, params)
    print("start_point (MG)", start_point)
    for i in range(int(integrate_time / step)):
        old_pointOOP = oop.copy()
        makeStep(start_point, dimension, diffFunc, params, step)
        oop = replace_var(start_point)
        ood = oop_to_ood(oop)

        if (old_pointOOP[2] > crossection and oop[2] < crossection and abs(oop[2]-old_pointOOP[2])<np.pi) or \
           (old_pointOOP[2] < crossection and oop[2] > crossection and abs(oop[2]-old_pointOOP[2])<np.pi):
            flag = 1
            dpt = replace_dpt(start_point, params)
            print("dpt", dpt)
            print("FIND CROSS")
            print("old", old_pointOOP)
            print("new", oop)
            makeStep(dpt, dimension_after_replace, diffPoincare, params, -(dpt[1] - crossection))
            print(i, start_point)
            dzeta_, phi_, teta_ = dpt
            omega = get_omega(dpt, params)
            # oop = np.array([omega * np.cos(dzeta_-phi_), omega * np.sin(dzeta_-phi_),phi_])
            ood = np.array([omega * np.cos(dzeta_-phi_), 
                            omega * np.sin(dzeta_-phi_),
                            dzeta_])
            break

        if np.isnan(start_point).any():
            print("NANANANAN_OG", i)
            break
        if np.isnan(oop).any():
            print("NANANANAN_OOP", i)
            break

        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = ood

    count_m[n] = i
    print("last point: ", start_point)
    # print(i, "time:", time.time()-start)
    if flag:
        return ood.copy()
    else:
        return np.array([np.nan,np.nan,np.nan])


figPoin = plt.figure(figsize=(6, 6))
axPoin = figPoin.add_subplot()
# axPoin.set_xlim(0,2*np.pi)
# axPoin.set_ylim(0,2*np.pi)

if __name__ == "__main__":
    start = time.time()
    phi = np.linspace(0, 2*np.pi, DISCR)
    points = []
    for _phi in phi:
        initial = np.array([0.01, _phi, 0.01])
        dzeta_, phi_, teta_ = initial
        omega = get_omega(initial, params)
        points.append(np.array([omega * np.cos(dzeta_-phi_), 
                                omega * np.sin(dzeta_-phi_), 
                                dzeta_]))

    start_points = np.array(points)
    arr = np.empty((len(phi), int((integrate_time / step) / skip_for_phase), dimension_after_replace))
    count_mas = np.empty(len(phi), dtype=int)

    mas_for_points = arr
    
    cross_points = []
    for number, i in enumerate(start_points):
        print(number, "start OOD",i)
        print("start DPT", [0.01,phi[number],0.01])
        cross_points.append(get_first_cross(i, mas_for_points[number], count_mas, number))

    cross_points = np.array(cross_points)
    print("колво точек на сечении", len(cross_points))
    cross_points = cross_points[~np.isnan(cross_points).any(axis=1)].T
    axPoin.plot(cross_points[0], cross_points[1], linestyle="", marker="o", markersize=0.5, color="black")

    plotter = pv.Plotter()

    line = pv.lines_from_points(start_points)
    plotter.add_mesh(line, color='black', line_width=3)

    colors = plt.cm.plasma(np.linspace(0, 1, len(mas_for_points)))
    points_number = 0
    for i, (traj, color) in enumerate(zip(mas_for_points, colors)):
        # if i == 0:
            # print(traj[-10:])
            # add_traj = traj[~np.isnan(traj).any(axis=1)]
            # print(mas_for_points.dtype)
            
            points_number += len(traj[:count_mas[i]])
            points = pv.PolyData(traj[:count_mas[i]])
            plotter.add_mesh(points, color=color[:3], point_size=1)

    bounds = plotter.bounds
    plotter.show_bounds(bounds=(bounds[0], bounds[1],   # X
                                bounds[2], bounds[3],   # Y
                                bounds[4], bounds[5]),     # Z
                        grid='back',        # сетка
                        location='outer',   # снаружи
                        all_edges=True      # все рёбра
    )
    print(points_number)
    print('result time', time.time()-start)
    plotter.show()

    plt.show()