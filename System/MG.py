'''
Тут рисуется просто фазовый портрет, 
интегрирование при этом в переменных M gamma
фазовый в переменных dzeta phi teta
'''
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


def main(start_point):
    mas_for_points = np.array([ [None] * dimension_after_replace for i in range(int((integrate_time / step) / skip_for_phase))])
    
    print("***********TIME for integrate:", (int(integrate_time / step)))

    start = time.time()
    # dpt = start_point
    # print(dpt)
    # for i in range(int(skip_time / step)):
    #     start_point = back_replace_var(dpt, params)
    #     makeStep(start_point, dimension, diffFunc, params, step)
    #     dpt = replace_var(start_point)
    #     print(dpt)
    #     norm_solution(dpt)
    #     if np.isnan(dpt[0]):
    #         print("NANANANAN")
    #         break

    # print("***********END warm time:", time.time() - start)

    start = time.time()
    En_mas = []
    geom_mas = []
    dpt = start_point
    for i in range(int(integrate_time / step)):
        # calc_energy_dpt(dpt, params)
        start_point = back_replace_var(dpt, params)
        # print("oldMG", start_point)
        # print("oldDPT", dpt)
        # dpt = replace_var(start_point)
        # print("newDPT", dpt)
        # start_point = back_replace_var(dpt, params)
        # print("newMG", start_point)
        # En_mas.append(calc_energy(start_point,params))
        # geom_mas.append(calc_geom_integral(start_point))
        # print("1",start_point)
        makeStep(start_point, dimension, diffFunc, params, step)
        # print("2",start_point)
        dpt = replace_var(start_point)
        # print(dpt)
        norm_solution(start_point)
        if np.isnan(start_point[0]):
            print("NANANANAN")
            break
        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = dpt
    end = time.time()
    print(mas_for_points[-10:])
    print("last point: ", start_point)
    print("***********TIME:", end-start)

    start = time.time()
    mas_for_points = mas_for_points.T
    mas_for_points = mas_for_points[mas_for_points != None]
    mas_for_points = mas_for_points.reshape((dimension_after_replace, int(len(mas_for_points) / dimension_after_replace)))

    # figEn, axEn = plt.subplots(1,2)
    # axEn[1].plot([i for i in range(len(En_mas))], En_mas)
    # axEn[0].plot([i for i in range(len(geom_mas))], geom_mas)

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