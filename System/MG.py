'''
Тут рисуется просто фазовый портрет, 
интегрирование при этом в переменных M gamma
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

def new_calc_energy(start_point, params):
    rho, mpar, I_1, I_2, a0, g, en = params
    I_3 = 0
    g_1, g_2, g_3 = start_point[2:]
    w_1, w_2 = start_point[:2]

    t10 = g_3 * g_3
    t20 = 1 - t10
    t40 = np.sqrt(t20)
    t50 = 1 / t40
    r1 = -rho*t50*g_1 - a0
    r2 = -rho*t50*g_2
    r3 = 0

    rg23 = r2*g_3 - r3*g_2
    rg13 = r1*g_3 - r3*g_1
    rg12 = r1*g_2 - r2*g_1
    sr23 = rg23*rg23
    sr13 = rg13*rg13
    sr12 = rg12*rg12

    I11 = I_1 + mpar *sr23
    I12 = -mpar*rg23*rg13
    I13 = mpar*rg23*rg12
    I21 = -mpar*rg23*rg13
    I22 = I_2 + mpar*sr13
    I23 = -mpar*rg13*rg12
    I31 = mpar*rg23*rg12
    I32 = -mpar*rg13*rg12
    I33 = I_3 + mpar*sr12

    # denom = I_1*I_2*I_3 + m * (I_1*I_2*sr12 + I_2*I_3*sr23 + I_1*I_3*sr13)

    # II11 = (I_2*I_3 + m*(I_2 * sr12 + I_3 * sr13)) / denom
    # II12 = I_3*m*rg13*rg23 / denom
    # II13 = -I_2*m*rg12*rg23 / denom
    # II21 = I_3*m*rg13*rg23 / denom
    # II22 = (I_1*I_3 + m*(I_1 * sr12 + I_3 * sr23)) / denom
    # II23 = I_1*m*rg12*rg13 / denom
    # II31 = -I_2*m*rg12*rg23 / denom
    # II32 = I_1*m*rg12*rg13 / denom
    # II33 = (I_1*I_2 + m*(I_1 * sr13 + I_2 * sr23)) / denom

    M1 = I11 * w_1 + I12 * w_2
    M2 = I21 * w_1 + I22 * w_2
	#M3 = I31 * w_1 + I32 * w_2;
    return (w_1 * M1 + w_2 * M2) / 2.0 + mpar*g*(rho*np.sqrt(1.0 - g_3 * g_3) + g_1 * a0)

def calc_energy_dpt(start_point,params, omega1, omega2):
    dzeta, phi, teta = start_point
    ksi = dzeta-phi
    ro, mpar, i1, i2, a0, g, en = params
    
    omega = np.sqrt(omega1 ** 2 + omega2 ** 2)

    E = ((omega ** 2) / 2) * (i1 * np.cos(ksi)**2 + i2 * np.sin(ksi)**2) + \
         (mpar / 2) * np.cos(teta)**2 * omega**2 * (ro * np.cos(phi+ksi) - \
                                                    -a0*np.sin(ksi)) ** 2 + \
         mpar * g * np.sin(teta) * (ro + a0 * np.sin(phi))
    print("E(dpt)", E)

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

    I = np.diag([i1, i2, 0])
    J = I + mpar * np.outer((np.cross(r, gamma)), (np.cross(r, gamma)))
    E = (1/2) * np.dot(omega, J @ omega) - mpar*g*np.dot(r,gamma)
    # E = (1/2) * mpar * np.dot(np.cross(r, omega), gamma) ** 2 + \
    #     (1/2) * np.dot(omega, I @ omega) - mpar * g * np.dot(r, gamma)
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
    start_point = back_replace_var(dpt, params)
    print("start_point (MG)",start_point)
    for i in range(int(integrate_time / step)):

        En_mas.append(calc_energy(start_point,params))
        geom_mas.append(calc_geom_integral(start_point))

        makeStep(start_point, dimension, diffFunc, params, step)
        dpt = replace_var(start_point)
        norm_solution(dpt)

        if np.isnan(start_point[0]):
            print("NANANANAN")
            break

        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = dpt
            
    end = time.time()
    # print(mas_for_points[-10:])
    print("last point: ", start_point)
    print("***********TIME:", end-start)

    start = time.time()
    mas_for_points = mas_for_points.T
    mas_for_points = mas_for_points[mas_for_points != None]
    mas_for_points = mas_for_points.reshape((dimension_after_replace, int(len(mas_for_points) / dimension_after_replace)))

    figEn, axEn = plt.subplots(1,2)
    axEn[0].set_ylim(0.9, 1.1)
    axEn[1].set_ylim(0, 2)
    axEn[1].plot([i for i in range(len(En_mas))], En_mas)
    axEn[0].plot([i for i in range(len(geom_mas))], geom_mas)

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