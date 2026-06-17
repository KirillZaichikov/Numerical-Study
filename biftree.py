'''
Код для построения однопараметрического дерева бифуркаций
'''


import numpy as np
import math as m
import matplotlib.pyplot as plt
from numba import jit


skip_time = 2000
len_time = 2000
DISCR = 10000 # disc param


params = np.array([0.00075, 0.56, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRM
border = 0.0009
min_border = 0.08


def main(params):
    fig, ax = plt.subplots(figsize=(12,8), dpi=400)
    ax.tick_params(axis='both', labelsize=18)
    # ax.set_xlim(params[0], border)
    mas_for_points = calc(params)
    mas_for_points = mas_for_points[~np.isnan(mas_for_points).any(axis=1)]
    mas_for_points = mas_for_points[~np.isinf(mas_for_points).any(axis=1)]
    mas_for_points = mas_for_points.T
    ax.plot(mas_for_points[0], mas_for_points[1], ',k', alpha=0.25)
    plt.savefig("biftree.png", dpi=400)
    plt.show()

@jit(nopython=True, cache=True)
def lerFRM_3D_map(state, res, params):
    eps, alpha, beta, p = params
    ksi0, eta0, teta0 = state

    # LOCAL 2
    ro = np.sqrt(ksi0**2+eta0**2)/p
    teta0 = teta0
    phi0 = np.atan2(eta0, ksi0) + teta0

    r = p * (ro / p) ** ((alpha+eps)/(alpha-eps))
    teta1 = ((beta + eps) / (alpha - eps)) * np.log(p/ro) + teta0
    phi1 = (((beta - eps) / (alpha - eps)) * np.log(p/ro) + phi0)

    while (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    while (phi1 < -m.pi):
        phi1 += 2 * m.pi
    ksi1 = r * p * np.cos(phi1-teta1)
    eta1 = r * p * np.sin(phi1-teta1)
    phi1 = phi1

    # ******** S3
    # ksi2 = ksi1 + eps * np.cos(2*phi1)
    # eta2 = eta1 + eps * np.sin(2*phi1)
    # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S2
    # ksi2 = ksi1 + eps * np.cos(phi1)
    # eta2 = eta1 + eps * np.sin(phi1)
    # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S1
    ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    while teta2>np.pi:
        teta2 -= 2 * np.pi
    while teta2<-np.pi:
        teta2 += 2 * np.pi

    return [ksi2, eta2, teta2]

@jit(nopython=True, cache=True)
def map(state, res, params):
    return lerFRM_3D_map(state,res,params)

@jit(nopython=True, cache=True)
def calc(params):
    # mas_for_points = np.zeros((DISCR*len_time, 2))
    # # mas_for_points = np.ndarray((DISCR*len_time, 2))
    # print('calc start:')
    # initial_point = [0.001,0.001,0]  # надо указывать
    # res = np.zeros(3)
    # step = (border - params[0]) / DISCR
    # new_params = params.copy()
    # for i in range(DISCR):
    #     new_params[0] += step
    #     print(new_params)

    #     for counter in range(skip_time):
    #         initial_point = map(initial_point, res, new_params)

    #     for counter in range(len_time):
    #         initial_point = map(initial_point, res, new_params)
    #         mas_for_points[counter+i*len_time] = [new_params[0], initial_point[0]]

    # mas_for_points = mas_for_points[~np.isnan(mas_for_points).any(axis=1)]
    # mas_for_points = mas_for_points[~np.isinf(mas_for_points).any(axis=1)]
    # print(mas_for_points)
    # mas_for_points = mas_for_points.T
    # ax.scatter(mas_for_points[0], mas_for_points[1], s=1, c="black")
    # ax.plot(mas_for_points[0], mas_for_points[1], linestyle='', marker='.', markersize=0.1, color="black")

    ####### ИДЕМ В ОБРАТНОМ НАПРАВЛЕНИИ От НАЧАЛЬНОЙ ТОЧКИ
    mas_for_points = np.zeros((DISCR*len_time, 2), dtype=np.float64)
    res = np.zeros(3)
    initial_point = [0.001,0.001,0]  # надо указывать
    params[0] = 0.1
    step = abs(min_border - params[0]) / DISCR
    for i in range(DISCR):
        params[0] -= step
        print(params)

        for counter in range(skip_time):
            initial_point = map(initial_point, res, params)

        for counter in range(len_time):
            initial_point = map(initial_point, res, params)
            mas_for_points[counter+i*len_time] = [params[0], initial_point[0]]

    # ax.plot(mas_for_points[0], mas_for_points[1], ',k', alpha=0.25)

    ####### ИДЕМ В ОБРАТНОМ НАПРАВЛЕНИИ
    # mas_for_points = np.ndarray((DISCR*len_time, 2))
    # # initial_point = [0.001,0.001,0]  # надо указывать
    # for i in range(DISCR):
    #     params[0] -= 0.2 / DISCR
    #     print(params)

    #     for counter in range(skip_time):
    #         initial_point = map(initial_point, res, params)

    #     for counter in range(len_time):
    #         initial_point = map(initial_point, res, params)
    #         mas_for_points[counter+i*len_time] = [params[0], initial_point[0]]

    # mas_for_points = mas_for_points[~np.isnan(mas_for_points).any(axis=1)]
    # mas_for_points = mas_for_points[~np.isinf(mas_for_points).any(axis=1)]
    # print(mas_for_points)
    # mas_for_points = mas_for_points.T
    # ax.scatter(mas_for_points[0], mas_for_points[1], s=0.1, c="skyblue")
    return mas_for_points

    

if __name__ == "__main__":
    main(params)