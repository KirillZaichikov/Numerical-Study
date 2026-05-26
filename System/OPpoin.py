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
# ax.set_ylim(0, 2*np.pi)
# ax.set_zlim(0, np.pi)

ax.set_xlabel("phi")
ax.set_ylabel("omega1")
ax.set_zlabel("omega2")

# @jit(nopython=True, cache=True)
def norm_solution(start_point):
    if start_point[2] > 2*np.pi:
        start_point[2] -= 2*np.pi
    elif start_point[2] < 0:
        start_point[2] += 2*np.pi

# @jit(nopython=True, cache=True)
def calc_energy_dpt(start_point,params):
    dzeta, phi, teta = start_point
    ksi = dzeta-phi
    ro, mpar, i1, i2, a0, g, en = params

    # omega = np.sqrt(2*(en - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
    #     (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
    #         mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2))
    
    
    E = ((omega ** 2) / 2) * (i1 * np.cos(ksi)**2 + i2 * np.sin(ksi)**2) + \
         (mpar / 2) * np.cos(teta)**2 * omega**2 * (ro * np.cos(phi+ksi) - \
                                                    -a0*np.sin(ksi)) ** 2 + \
         mpar * g * np.sin(teta) * (ro + a0 * np.sin(phi))
    # print("E(dpt)", E)

# @jit(nopython=True, cache=True)
def back_replace_var_dpt(start_point, params):
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

# @jit(nopython=True, cache=True)
def replace_var(start_point):
    '''
    MG -> omega1, omega2, phi
    '''
    omega1, omega2, _g1, _g2, _g3 = start_point
    _phi = np.atan2(_g1, _g2)
    point = np.array([omega1, omega2, _phi])
    norm_solution(point)
    return point.copy()

# @jit(nopython=True, cache=True)
def back_replace_var(start_point, params):
    '''
    omega1,omega2,phi -> OmegaGamma
    '''
    omega1, omega2, phi = start_point
    ro, mpar, i1, i2, a0, g, en = params

    ksi = np.atan2(omega2, omega1)
    dzeta = ksi + phi
    omega = np.sqrt(omega1 ** 2 + omega2 ** 2)

    # Результат зависит от начального приближения
    # teta0 = np.pi-eps
    teta0 = 0 + 0.1
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

# @jit(nopython=True, cache=True)
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
    # print("Energy", E)
    return E

# @jit(nopython=True, cache=True)
def calc_geom_integral(start_point):
    geom = np.linalg.norm(start_point[2:])
    # print("geom", geom)
    return geom

# @jit(nopython=True, cache=True)
def get_omega(start_point, params):
    ro, mpar, i1, i2, a0, g, eps = params
    dzeta, phi, teta = start_point

    omega2 = 2*(eps - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
        (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
            mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2)
    omega = m.sqrt(omega2)
    if omega2 < 0: print("ERROR: omega less than 0")
    return omega

# @jit(nopython=True, cache=True)
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

# @jit(nopython=True, cache=True)
def main(start_point):
    mas_for_points = np.array([ [None] * dimension_after_replace for i in range(int((integrate_time / step) / skip_for_phase))])
    mas_for_cross = []
    
    print("***********TIME for integrate:", (int(integrate_time / step)))

    start = time.time()
    En_mas = []
    geom_mas = []
    oop_start = start_point.copy()
    oop = start_point
    start_point = back_replace_var(oop, params)
    print("start_point (MG)", start_point)
    # omega = get_omega(oop, params)
    # gam = np.array([0,0.01, 1-0.01])
    # nor = np.linalg.norm(gam)
    # for i in range(len(gam)):
    #     gam[i] = gam[i] / nor
    # print(np.linalg.norm(gam))
    # start_point = np.array([0,0,1,omega * np.cos(dpt[0]-dpt[1]), omega * np.sin(dpt[0]-dpt[1])])
    # start_point = np.array([omega * np.cos(dpt[0]-dpt[1]), omega * np.sin(dpt[0]-dpt[1]), gam[0], gam[1], gam[2]])
    for i in range(int(integrate_time / step)):

        old_point = start_point.copy()
        old_pointOOP = oop.copy()
        makeStep(start_point, dimension, diffFunc, params, step)

        En_mas.append(calc_energy(start_point,params))
        geom_mas.append(calc_geom_integral(start_point))
        oop = replace_var(start_point)

        if (old_pointOOP[2] > crossection and oop[2] < crossection and abs(oop[2]-old_pointOOP[2])<np.pi) or \
           (old_pointOOP[2] < crossection and oop[2] > crossection and abs(oop[2]-old_pointOOP[2])<np.pi):
            flag = 1
            dpt = replace_dpt(start_point, params)
            print("FIND CROSS")
            print("dpt", dpt)
            print("old", old_pointOOP)
            print("new", oop)
            makeStep(dpt, dimension_after_replace, diffPoincare, params, -(dpt[1] - crossection))
            print(i, dpt)
            dzeta_, phi_, teta_ = dpt
            omega = get_omega(dpt, params)
            oop = np.array([omega * np.cos(dzeta_-phi_), omega * np.sin(dzeta_-phi_),phi_])
            mas_for_cross.append(oop)
            start_point = back_replace_var_dpt(dpt, params)
            print("cross oop", oop)

        if np.isnan(start_point[0]):
            print("NANANANAN")
            break

        if i % skip_for_phase == 0:
            mas_for_points[i // skip_for_phase] = oop

    end = time.time()
    # print(mas_for_points[-10:])
    print("old_pointOOP:", old_pointOOP)
    print("old_point:", old_point)
    print("last point: ", start_point)
    print("***********TIME:", end-start)

    start = time.time()
    mas_for_cross = np.array(mas_for_cross).T
    mas_for_points = mas_for_points.T
    # mas_for_points = mas_for_points[mas_for_points != None]
    # mas_for_points = mas_for_points.reshape((dimension_after_replace, int(len(mas_for_points) / dimension_after_replace)))
    return mas_for_cross, mas_for_points, En_mas, geom_mas, oop_start, old_pointOOP
    

if __name__ == "__main__":
    initial_point = np.array([0.01,0.06346651825433926,0.01])
    dzeta, phi, teta = initial_point
    print("start_point (dpt)", initial_point)

    omega = get_omega(initial_point, params)
    new_initial = np.array([omega * np.cos(dzeta-phi), omega * np.sin(dzeta-phi), phi])

    start_point = np.array([new_initial])
    print("start_point (OOP)", start_point)
    for i in start_point:
        mas_for_cross, mas_for_points, En_mas, geom_mas, oop_start, old_pointOOP = main(i)

    figEn, axEn = plt.subplots(1,2)
    axEn[0].set_ylim(0.9, 1.1)
    axEn[1].set_ylim(0, 2)
    axEn[1].plot([i for i in range(len(En_mas))], En_mas)
    axEn[0].plot([i for i in range(len(geom_mas))], geom_mas)

    ax.plot(mas_for_points[2], mas_for_points[0], mas_for_points[1], linestyle="", marker="o", markersize=0.1, color="black")
    ax.scatter(oop_start[2], oop_start[0], oop_start[1], c="red")
    ax.scatter(old_pointOOP[2], old_pointOOP[0], old_pointOOP[1], c="blue")
    ax.scatter(mas_for_cross[2],mas_for_cross[0],mas_for_cross[1], c="magenta")

    figPoin, axPoin = plt.subplots()
    axPoin.plot(mas_for_cross[0], mas_for_cross[1], linestyle="", marker="o", markersize=2, color="black")

    print("*********** calc is ending")
    # plt.savefig("plot.pdf", dpi=300)
    plt.show()