'''
Вычисляет ляпуновские показатели в координатах кузнецова
'''

import numpy as np
import matplotlib.pyplot as plt
import math as m
from utils.linal import *
import tqdm

omega_y = []

def calc_r_vec(a1, a2, h, gamma, d):
    gamma1, gamma2, gamma3 = gamma[0], gamma[1], gamma[2]
    r1 = -(a1 * gamma1) / gamma3 # КУЗНЕЦОВ
    r2 = -(a2 * gamma2) / gamma3 # КУЗНЕЦОВ
    r3 = -h + (1/2)*((a1 * gamma1 ** 2 + a2 * gamma2 ** 2)/(gamma3 ** 2)) # КУЗНЕЦОВ
    r = np.array([r1,r2,r3])
    return r

def calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M):
    Q = np.array([[np.cos(d), np.sin(d), 0],
                 [-np.sin(d),np.cos(d),0],
                 [0,0,1]])
    Q = np.array([[np.cos(d), np.sin(d), 0],
                 [-np.sin(d),np.cos(d),0],
                 [0,0,1]])
    # Формула бизяева в maple
    A = Q @ np.diag([I1,I2,I3]) @ Q.T + np.identity(3) * np.dot(r,r) - np.outer(r, r)
    A_inv = np.linalg.inv(A)
    omega = A_inv @ M
    return omega

def sys(state, res, params, H): # масса 1
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    M = np.array([state[0], state[1], state[2]])
    gamma = np.array([state[3], state[4], state[5]])

    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    # print(r, omega)

    # tmp = M.copy()
    # for k in range(3):
    #     M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))

    # tmp = gamma.copy()
    # for k in range(3):
    #     gamma[k] = gamma[k]/(np.sqrt(np.dot(tmp, tmp)))


    dgamma = np.cross(gamma, omega)
    Jr = np.array([[-a1/gamma[2], 0, a1 * gamma[0]/(gamma[2]**2)],
                   [0, -a2/gamma[2], a2 * gamma[1]/(gamma[2]**2)],
                   [a1*gamma[0]/(gamma[2]**2), a2*gamma[1]/(gamma[2]**2), -(a1*gamma[0]**2 + a2 * gamma[1]**2)/(gamma[2]**3)]])
    dr = Jr @ dgamma

    dM = np.cross(M, omega) + np.cross(dr, np.cross(omega, r)) + g0 * np.cross(r, gamma)

    res[0], res[1], res[2], res[3], res[4], res[5] = dM[0], dM[1], dM[2], dgamma[0], dgamma[1], dgamma[2]
    # res[0], res[1], res[2], res[3], res[4], res[5] = m1_, m2_, m3_, gamma1_, gamma2_, gamma3_

def sysSlave(state, res, params, H): # масса 1
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    M = np.array([state[0], state[1], state[2]])
    gamma = np.array([state[3], state[4], state[5]])

    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    # print(r, omega)

    ###################экспериментально привожу к первым интегралам
    # tmp = M.copy()
    # for k in range(3):
    #     M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))

    # tmp = gamma.copy()
    # for k in range(3):
    #     gamma[k] = gamma[k]/(np.sqrt(np.dot(tmp, tmp)))

    # En = (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma)
    # print("energy", En)
    # print("geom", np.dot(gamma, gamma))

    dgamma = np.cross(gamma, omega)
    Jr = np.array([[-a1/gamma[2], 0, a1 * gamma[0]/(gamma[2]**2)],
                   [0, -a2/gamma[2], a2 * gamma[1]/(gamma[2]**2)],
                   [a1*gamma[0]/(gamma[2]**2), a2*gamma[1]/(gamma[2]**2), -(a1*gamma[0]**2 + a2 * gamma[1]**2)/(gamma[2]**3)]])
    dr = Jr @ dgamma

    dM = np.cross(M, omega) + np.cross(dr, np.cross(omega, r)) + g0 * np.cross(r, gamma)

    res[0], res[1], res[2], res[3], res[4], res[5] = dM[0], dM[1], dM[2], dgamma[0], dgamma[1], dgamma[2]
    # res[0], res[1], res[2], res[3], res[4], res[5] = m1_, m2_, m3_, gamma1_, gamma2_, gamma3_

def dverkStep(val, dimension, diffFunc, params, step, H=1) -> None:
    k1, k2, k3, k4, k5, k6, k7, k8 = np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), \
        np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension)
    arg = np.zeros(dimension)
    diffFunc(val, k1, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step / 6. * k1[j]
    
    diffFunc(arg, k2, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j])
    
    diffFunc(arg, k3, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j])
    
    diffFunc(arg, k4, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j])
    
    diffFunc(arg, k5, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j])
    
    diffFunc(arg, k6, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j])
    
    diffFunc(arg, k7, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j])
    
    diffFunc(arg, k8, params, H)
    for j in range(dimension):
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j])


if __name__ == "__main__":
    x1, x2 = [], []

    # Начальная точка из хаоса
    M = np.array([176.39671645238903, -168.85257112196294, -94.877918535603939])
    gamma = np.array([0.69524430378189617, 0.40429993079272591, -0.59428864887075838])
    gamma = -gamma

    # Параметры для вычислений
    time_skip = 10
    time_skip_clv = 10
    time = 10
    step = 0.01
    lyap_num = 6
    dimension = 6
    eps = 0.0001

    # Параметры из хаоса
    params = [0.485, 2, 6, 7, 9, 4, 1, 752, 100] # d I1 I2 I3 a1 a2 h E g0
    # params = [0.2, 5, 6, 7, 9, 4, 1, 555, 100]
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    # Замена начальных условий из хаоса на значения кузнецова
    QC = np.array([[np.cos(d), np.sin(d), 0],
                 [-np.sin(d),np.cos(d),0],
                 [0,0,1]])
    M = QC @ M
    gamma = QC @ gamma

    # посчитаем уровень энергии
    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    En = (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma)
    print(En)

    new_norm = np.zeros(dimension)
    P = [0 for i in range(dimension)]
    slaveTrajectory = np.identity(dimension) * eps

    main_traj = np.array([*M, *gamma])
    for i in range (int(time_skip/step)):
        dverkStep(main_traj, dimension, sys, params, step)

    print(main_traj)

    for i in range (int(time_skip_clv/step)):
        for j in range(lyap_num):
            slaveTrajectory[j] += main_traj

        dverkStep(main_traj, dimension, sys, params, step)

        for j in range(lyap_num):
            dverkStep(slaveTrajectory[j], dimension, sys, params, step)
            slaveTrajectory[j] -= main_traj
        ortVecs(slaveTrajectory, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(slaveTrajectory[j])
            slaveTrajectory[j] = (slaveTrajectory[j] / new_norm[j]) * eps
        # for j in range(lyap_num):
        #     print(np.linalg.norm(slaveTrajectory[j]))

    print(main_traj)
    for row in slaveTrajectory:
        # print(np.linalg.norm(row))
        print(" ".join(f"{val:10.8f}" for val in row))

    for i in range(int(time/step)):
        for j in range(lyap_num):
            slaveTrajectory[j] += main_traj

        dverkStep(main_traj, dimension, sys, params, step)

        for j in range(lyap_num):
            dverkStep(slaveTrajectory[j], dimension, sys, params, step)
            slaveTrajectory[j] -= main_traj
            # print("1", slaveTrajectory[j])
            # print("2", slaveTrajectory[j])
            
        ortVecs(slaveTrajectory, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(slaveTrajectory[j])
            slaveTrajectory[j] = (slaveTrajectory[j] / new_norm[j]) * eps
            P[j] += m.log(new_norm[j] / eps)

    print(main_traj)

    for i in range(lyap_num):
        print(f"L{i+1}: ", P[i] / time)

    # r = calc_r_vec(a1, a2, h, gamma, d)
    # omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    # En = (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma)
    # print(En)
    plt.show()