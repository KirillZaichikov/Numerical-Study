'''
фазовый портрет и интегралы в координатах кузнецова
'''

import numpy as np
import matplotlib.pyplot as plt

# Начальная точка из хаоса
# M = np.array([-34.697048094358898, -62.831273385025732, 62.2309106274670827])
# gamma = np.array([-0.16059001362453112, -0.29080499936213250, -0.94320904356884894])
# M = np.array([176.39671645238903, -168.85257112196294, -94.877918535603939])
# gamma = np.array([0.69524430378189617, 0.40429993079272591, -0.59428864887075838])
# это я протянул почти к рождению УЖЕ ПОВЕРНУТЫЕ!!!
M = np.array([-59.79105323, -40.95557389,  61.6906867]) 
gamma = np.array([0.27751469, 0.19009154, 0.94172756])
# gamma = -gamma

time_skip = 0
time = 100
step = 0.0025
# Параметры из хаоса
#         d     I1 I2 I3 a1 a2  h   E    g0
# params = [0.489, 2, 6, 7, 9, 4, 1, 752, 100] 
params = [0.433, 2, 6, 7, 9, 4, 1, 740, 100] 
# params = [0.2  , 5, 6, 7, 9, 4, 1, 555, 100]
d = params[0]
I1, I2, I3 = params[1], params[2], params[3]
a1, a2 = params[4], params[5]
h, E, g0 = params[6], params[7], params[8]


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
    # print("JQ_inv\n", A_inv)
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

    ############# ЭКСПЕРИМЕНТАЛЬНО
    # try:
    #     tmp = M.copy()
    #     for k in range(3):
    #         M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))
    # except:
    #     pass

    # tmp = gamma.copy()
    # for k in range(3):
    #     gamma[k] = gamma[k]/(np.sqrt(np.dot(tmp, tmp)))

    #####################

    # print('r', r, "\n", 'omega', omega)

    
    dgamma = np.cross(gamma, omega)
    # print("dgamma", dgamma)
    Jr = np.array([[-a1/gamma[2], 0, a1 * gamma[0]/(gamma[2]**2)],
                   [0, -a2/gamma[2], a2 * gamma[1]/(gamma[2]**2)],
                   [a1*gamma[0]/(gamma[2]**2), a2*gamma[1]/(gamma[2]**2), -(a1*gamma[0]**2 + a2 * gamma[1]**2)/(gamma[2]**3)]])
    dr = Jr @ dgamma
    # print("dr", dr)
    dM = np.cross(M, omega) + np.cross(dr, np.cross(omega, r)) + g0 * np.cross(r, gamma)
    # print("dM", dM)
    res[0], res[1], res[2], res[3], res[4], res[5] = dM[0], dM[1], dM[2], dgamma[0], dgamma[1], dgamma[2]
    # res[0], res[1], res[2], res[3], res[4], res[5] = m1_, m2_, m3_, gamma1_, gamma2_, gamma3_

def dverkStep(val, dimension, diffFunc, params, step, H=1) -> None:
    k1, k2, k3, k4, k5, k6, k7, k8 = np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), \
        np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension)
    arg = np.zeros(dimension)
    diffFunc(val, k1, params, H)
    # print(val)
    # print(k1)
    for j in range(dimension):
        arg[j] = val[j] + step / 6. * k1[j]
        # print(arg[j])

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

def calcIntegrals(M, g0, gamma):
    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    print("(calcIntegrals) r", r, "omega", omega)
    return (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma), np.dot(gamma, gamma)

if __name__ == "__main__":
    x1, x2 = [], []
    # Замена начальных условий из хаоса на значения кузнецова (ЕСЛИ НАДО!!!)
    # QK = np.array([[np.cos(-d), np.sin(-d), 0],
    #              [-np.sin(-d),np.cos(-d),0],
    #              [0,0,1]])
    # QC = np.array([[np.cos(d), np.sin(d), 0],
    #              [-np.sin(d),np.cos(d),0],
    #              [0,0,1]])
    # M = QC @ M
    # gamma = QC @ gamma

    # Посчитаем r и omega для приведения M к необходимому уровню энергии (ЕСЛИ НАДО!!!)
    # r = calc_r_vec(a1, a2, h, gamma, d)
    # # print('вектор r ', r)
    # omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    # # print('вектор omega ', omega)
    # tmp = M.copy()
    # print(tmp)
    # for k in range(3):
    #     M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))
    # print('M после нормировки', M)
    
    # Для новых m и gamma посчитаем новые r и omega
    # r = calc_r_vec(a1, a2, h, gamma, d)
    # omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    # En = (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma)
    print("Энергия после поворота и приведения", *calcIntegrals(M, g0, gamma))

    En_list = []
    Geom_list = []
    M1_list, M2_list, M3_list, gamma1_list, gamma2_list, gamma3_list = [], [], [], [], [], []
    old_M = M.copy()
    main_traj = np.array([*M, *gamma])

    for i in range (int(time_skip/step)):
        dverkStep(main_traj, 6, sys, params, step)
    M[0], M[1], M[2], gamma[0], gamma[1], gamma[2] = main_traj[0], main_traj[1], main_traj[2], main_traj[3], main_traj[4], main_traj[5]

    for i in range (int(time/step)):

        main_traj = np.array([*M, *gamma])
        dverkStep(main_traj, 6, sys, params, step)
        M[0], M[1], M[2], gamma[0], gamma[1], gamma[2] = main_traj[0], main_traj[1], main_traj[2], main_traj[3], main_traj[4], main_traj[5]
        # print(main_traj)

        r = calc_r_vec(a1, a2, h, gamma, d)
        omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
        En = (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma)
        
        x1.append(i)
        En_list.append(En)
        Geom_list.append(np.dot(gamma, gamma))
        M1_list.append(M[0])
        M2_list.append(M[1])
        M3_list.append(M[2])
        gamma1_list.append(gamma[0])
        gamma2_list.append(gamma[1])
        gamma3_list.append(gamma[2])
        # M = QK@M
        # gamma = QK @ gamma

    print(En)
    fig, ax = plt.subplots(2,2)
    ax[0][0].set_title("Energy")
    ax[0][0].plot(x1, En_list)
    ax[0][0].set_ylim(En-10, En+10)

    ax[0][1].set_title("Geom Integral")
    ax[0][1].plot(x1, Geom_list)
    ax[0][1].set_ylim(0.999, 1.001)

    ax[1][0].set_title("M1")
    ax[1][0].plot(x1, M1_list)

    ax[1][1].set_title("gamma1")
    ax[1][1].plot(x1, gamma1_list)

    fig_phase = plt.figure()
    ax1 = fig_phase.add_subplot(2, 2, 1, projection='3d')
    ax1.scatter(0,0,0, c='red', s=5)
    ax1.plot(M1_list, M2_list, gamma1_list)

    ax2 = fig_phase.add_subplot(2, 2, 2, projection='3d')
    ax2.scatter(0,0,0, c='red', s=5)
    ax2.plot(M1_list, M2_list, M3_list)

    ax3 = fig_phase.add_subplot(2, 2, 3, projection='3d')
    ax3.scatter(0,0,0, c='red', s=5)
    ax3.plot(M2_list, gamma1_list, gamma2_list)

    ax4 = fig_phase.add_subplot(2, 2, 4, projection='3d')
    ax4.scatter(0,0,0, c='red', s=5)
    ax4.plot(gamma3_list, gamma1_list, gamma2_list)

    plt.show()
