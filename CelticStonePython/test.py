import numpy as np
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    x1, x2 = [], []

    # Начальная точка из хаоса
    # M = np.array([-34.697048094358898, -62.831273385025732, 62.2309106274670827])
    # gamma = np.array([-0.16059001362453112, -0.29080499936213250, -0.94320904356884894])
    M = np.array([-0.6033344705e-2, 0.3002787660e-1, 0.1866956664])
    gamma = np.array([-8.731875373*10^(-7), -0.1247953839e-5, -0.414878368e-5])

    time_skip = 50
    time = 10
    step = 0.0025
    # Параметры из хаоса
    #         d     I1 I2 I3 a1 a2  h   E    g0
    params = [0.480, 2, 6, 7, 9, 4, 1, 780, 100] 
    # params = [0.2  , 5, 6, 7, 9, 4, 1, 555, 100]
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    # Замена начальных условий из хаоса на значения кузнецова
    QK = np.array([[np.cos(-d), np.sin(-d), 0],
                 [-np.sin(-d),np.cos(-d),0],
                 [0,0,1]])
    QC = np.array([[np.cos(d), np.sin(d), 0],
                 [-np.sin(d),np.cos(d),0],
                 [0,0,1]])
    M = QC @ M
    # print('M после поворота ', M)
    gamma = QC @ gamma
    # print('gamma после поворота', gamma)

    # Посчитаем r и omega для приведения M к необходимому уровню энергии
    r = calc_r_vec(a1, a2, h, gamma, d)
    # print('вектор r ', r)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    # print('вектор omega ', omega)
    tmp = M.copy()
    print(tmp)
    for k in range(3):
        M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))
    print('M после нормировки', M)
    
    # Для новых m и gamma посчитаем новые r и omega
    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    En = (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma)
    print("Энергия после поворота и приведения",En)
    print(r, omega)

    En_list = []
    Geom_list = []
    M1_list, M2_list, gamma1_list = [], [], []
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
        gamma1_list.append(gamma[0])
        # M = QK@M
        # gamma = QK @ gamma

    print(En)
    fig, ax = plt.subplots(2,2)
    ax[0][0].set_title("Energy")
    ax[0][0].plot(x1, En_list)
    # ax[0][0].set_ylim(460, 520)
    # ax[0][0].set_ylim(712, 772)
    # ax[0][0].set_ylim(1350, 1420)
    ax[0][1].set_title("Geom Integral")
    ax[0][1].plot(x1, Geom_list)
    ax[0][1].set_ylim(0.999, 1.001)
    ax[1][0].set_title("M1")
    ax[1][0].plot(x1, M1_list)
    ax[1][1].set_title("gamma1")
    ax[1][1].plot(x1, gamma1_list)

    fig_phase = plt.figure()
    ax_phase = fig_phase.add_subplot(projection='3d')
    ax_phase.plot(M1_list, M2_list, gamma1_list)
    plt.show()
