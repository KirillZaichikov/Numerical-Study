'''
Фазовый портрет для кельтского камня в координатах хаоса
'''

import numpy as np
import matplotlib.pyplot as plt

omega_y = []
# def sys(res, m1, m2, m3, gamma1, gamma2, gamma3, d, I1, I2, I3, a1, a2, h, E, g0): # масса 1, gamma3 считаем из геом интеграла


# ПОКА НЕ НАДО
def calc_gamma3(gamma1, gamma2):
    # из геом интеграла g1^2 + g2^2 + g3^2 = 1 выразим gamma3 пусть с плюсом
    # print("проверка геометрического интеграла", gamma1**2+gamma2**2+m.sqrt(1-gamma1**2-gamma2**2)**2)
    return -m.sqrt(1-gamma1**2-gamma2**2)

# ПОКА НЕ НАДО
def calc_m3_and_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M):
    m1, m2 = M[0], M[1]
    r1,r2,r3 = r[0], r[1], r[2]

    A = np.diag([I1,I2,I3]) + np.identity(3) * np.dot(r,r) - np.outer(r, r)
    A_inv = np.linalg.inv(A)
    a11, a12, a13 = A_inv[0][0], A_inv[0][1], A_inv[0][2]
    a21, a22, a23 = A_inv[1][0], A_inv[1][1], A_inv[1][2]
    a31, a32, a33 = A_inv[2][0], A_inv[2][1], A_inv[2][2]
    # print('gamma', gamma)
    B = a13 * m1 + a31 * m1 + a23 * m2 + a32 * m2
    C = a11 * m1 ** 2 + (a12+a21)*m1*m2 + a22 * m2 ** 2
    D = (B/2)**2-4*(a33/2)*(C/2-g0 * np.dot(r, gamma) - E)
    # print(D)
    m3 = (-B/2 + (D) ** (1/2)) / a33                                              # взял m3 тот что с плюсом
    M = np.array([m1,m2,m3])
    omega = A_inv @ np.array([m1, m2, m3])
    # print('omega', omega)
    omega1 = omega[0]
    omega2 = omega[1]
    omega3 = omega[2]
    # omega_y.append(omega1)
    #проверка
    print("(calc_m3_func) OMEGA", omega)
    return m3, omega1, omega2, omega3
    # pass

# ВЕКТОР R СХОДИТСЯ
def calc_r_vec(a1, a2, h, gamma, d):
    # ЭТО ТО ЧТО Я ВБИВАЛ РУКАМИ ОНО ТОЖЕ РАБОТАЕТ
    # gamma1, gamma2, gamma3 = gamma[0], gamma[1], gamma[2]
    # r1 = (- 1 / gamma3) * \
    #     ((a1*np.cos(d)**2 + a2*np.sin(d)**2)*gamma1 + \
    #      (a1*np.cos(d)*np.sin(d)-a2*np.cos(d)*np.sin(d))*gamma2)
    # r2 = (- 1 / gamma[2]) * \
    #     ((a1*np.cos(d)*np.sin(d)-a2*np.cos(d)*np.sin(d))*gamma1 + \
    #      (a1*np.sin(d)**2+a2*np.cos(d)**2)*gamma2)
    # r3 = -h + (1/(2*gamma3**2)) * \
    #     (((a1*np.cos(d)**2 + a2*np.sin(d)**2)*gamma1 + \
    #       (a1 * np.cos(d)*np.sin(d)-a2*np.cos(d)*np.sin(d))*gamma2)*gamma1 + \
    #     ((a1*np.cos(d)*np.sin(d)-a2*np.cos(d)*np.sin(d))*gamma1+\
    #      (a1*np.sin(d)**2+a2*np.cos(d)**2)*gamma2)*gamma2)
    # r = np.array([r1,r2,r3])

    # ФОРМУЛА ИЗ МЭПЛА
    r = np.ndarray(3)
    a = [a1, a2]
    r[0] = -0.1e1 / gamma[2] * ((a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1])
    r[1] = -0.1e1 / gamma[2] * ((a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1])
    r[2] = -h + (((a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1]) * gamma[0] + ((a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1]) * gamma[1]) * pow(gamma[2], -0.2e1) / 0.2e1
    return r

# Вектор omega СХОДИТСЯ
def calc_omega(I1, I2, I3, r, gamma, M):
    r1,r2,r3 = r[0], r[1], r[2]
    # Формула бизяева в maple
    A = np.diag([I1,I2,I3]) + np.identity(3) * np.dot(r,r) - np.outer(r, r)
    A_inv = np.linalg.inv(A)
    omega = A_inv @ M
    # omega1 = A_inv[0][0] * m1 + A_inv[0][1] * m2 + A_inv[0][2] * m3
    # omega2 = A_inv[1][0] * m1 + A_inv[1][1] * m2 + A_inv[1][2] * m3
    # omega3 = A_inv[2][0] * m1 + A_inv[2][1] * m2 + A_inv[2][2] * m3
    # omega = np.array([omega1, omega2, omega3])
    print('(calc_omega)OMEGA', omega)
    # omega1 = omega[0]
    # omega2 = omega[1]
    # omega3 = omega[2]
    #проверка
    # print("проверка интеграла энергии", 1/2 * np.dot(M, omega) - g0 * np.dot(r, gamma)-E)
    # return omega1, omega2, omega3
    return omega

def sys(state, res, params, omega, r, H): # масса 1
    m1, m2, m3 = state[0], state[1], state[2]
    gamma1, gamma2, gamma3 = state[3], state[4], state[5]
    r1, r2, r3 = r[0], r[1], r[2]

    g0 = params[8]
    a = [params[4], params[5]]
    # a1, a2 = params[4], params[5]  

    M = np.array([state[0], state[1], state[2]])
    gamma = np.array([state[3], state[4], state[5]])

    # Матрица соотношения между M и omega
    # A = np.array([[I1 * (m.cos(d))**2 + I2*(m.sin(d))**2 + r2 ** 2 + r3 ** 2, 
    #               (I1-I2) * m.cos(d) * m.sin(d) - r1 * r2,
    #               - r1 * r2],
    #              [(I1-I2) * m.cos(d) * m.sin(d) - r1 * r2,
    #               I1 * (m.sin(d))**2 + I2*(m.cos(d))**2 + r1 ** 2 + r3 ** 2,
    #               - r2 * r3],
    #              [- r1 * r3, - r2 * r3, I3 + r1 ** 2 + r2 ** 2]])
    # A_inv = np.linalg.inv(A)
    # print('проверка обратной матрицы', A @ A_inv)
    # a11, a12, a13 = A_inv[0][0], A_inv[0][1], A_inv[0][2]
    # a21, a22, a23 = A_inv[1][0], A_inv[1][1], A_inv[1][2]
    # a31, a32, a33 = A_inv[2][0], A_inv[2][1], A_inv[2][2]

    # B = a13 * m1 + a31 * m1 + a23 * m2 + a32 * m2
    # C = a11 * m1 ** 2 + (a12+a21)*m1*m2 + a22 * m2 ** 2
    # D = (B/2)**2-4*a33/2*(C/2-g0 * np.dot(r, gamma) - E)
    # m3 = (-B/2 + m.sqrt(D)) / a33                                              # взял m3 тот что с плюсом
    # M = np.array([m1,m2,m3])
    # omega = A_inv @ np.array([m1, m2, m3])
    # omega1 = omega[0]
    # omega2 = omega[1]
    # omega3 = omega[2]
    # omega = np.array([omega1,omega2,omega3])
    # omega_y.append(omega1)
    # #проверка
    # print("проверка интеграла энергии", 1/2 * np.dot(M, omega) - g0 * np.dot(r, gamma)-E)

    # # из интеграла энергии выразим m3
    # m3 = (2 * E * 2 * g0 * (r1 * gamma1 + r2 * gamma2 + r3 * gamma3) - m1 * omega1 + m2 * omega2) / omega3

    #ТАК ПИСАЛ Я

    # Вычисление p и omegar 
    # r1_ = (- a1 * ((gamma2 * omega3 - gamma3 * omega2) * gamma3 - 
    #               (gamma1 * omega2 - gamma2 * omega1) * gamma1)) / (gamma3**2) 
    # r2_ = (- a2 * ((gamma3 * omega1 - gamma1 * omega3) * gamma3 - 
    #               (gamma1 * omega2 - gamma2 * omega1) * gamma2)) / (gamma3**2) 
    # r3_ = (a1 * gamma1 * (gamma2 * omega3 - gamma3 * omega2) + 
    #        a2 * gamma2 * (gamma3 * omega1 - gamma1 * omega3)) / (gamma3**2) - \
    #       (a1 * gamma1 ** 2 + a2 * gamma2 ** 2)*(gamma1 * omega2 - gamma2 * omega1)/(gamma3 ** 3)
    # p = r1 * r1_ + r2 * r2_ + r2 * r3_
    # omegar = r1_ * omega1 + r2_ * omega2 + r3_ * omega3

    # m1_ = m2 * omega3 - m3 * omega2 + omega1 * p - r1 * omegar + g0 * (r2 * gamma3 - gamma2 * r3)
    # m2_ = m3 * omega1 - m1 * omega3 + omega2 * p - r2 * omegar + g0 * (r3 * gamma1 - gamma3 * r1)
    # m3_ = m1 * omega2 - m2 * omega1 + omega3 * p - r3 * omegar + g0 * (r1 * gamma2 - gamma1 * r2)
    # gamma1_ = gamma2 * omega3 - gamma3 * omega2
    # gamma2_ = gamma3 * omega1 - gamma1 * omega3
    # gamma3_ = gamma1 * omega2 - gamma2 * omega1

    # А ТАК БИЗЯЕВ
    dgamma = np.cross(gamma, omega)

    dr = np.ndarray(3)
    dr[0] = -0.1e1 / gamma[2] * (a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * dgamma[0] - 0.1e1 / gamma[2] * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * dgamma[1] + pow(gamma[2], -0.2e1) * ((a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1]) * dgamma[2]
    dr[1] = -0.1e1 / gamma[2] * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * dgamma[0] - 0.1e1 / gamma[2] * (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * dgamma[1] + pow(gamma[2], -0.2e1) * ((a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1]) * dgamma[2]
    dr[2] = (0.2e1 * (a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + 0.2e1 * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1]) * pow(gamma[2], -0.2e1) * dgamma[0] / 0.2e1 + (0.2e1 * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + 0.2e1 * (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1]) * pow(gamma[2], -0.2e1) * dgamma[1] / 0.2e1 - (((a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1]) * gamma[0] + ((a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1]) * gamma[1]) * pow(gamma[2], -0.3e1) * dgamma[2]

    dM = np.cross(M, omega) + np.cross(dr, np.cross(omega, r)) + g0 * np.cross(r, gamma)

    res[0], res[1], res[2], res[3], res[4], res[5] = dM[0], dM[1], dM[2], dgamma[0], dgamma[1], dgamma[2]
    # res[0], res[1], res[2], res[3], res[4], res[5] = m1_, m2_, m3_, gamma1_, gamma2_, gamma3_

def dverkStep(val, dimension, diffFunc, params, step, omega, r, H=1) -> None:
    k1, k2, k3, k4, k5, k6, k7, k8 = np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), \
        np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension)
    arg = np.zeros(dimension)
    diffFunc(val, k1, params, omega, r, H)
    for j in range(dimension):
        arg[j] = val[j] + step / 6. * k1[j]
    
    diffFunc(arg, k2, params, omega, r, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j])
    
    diffFunc(arg, k3, params, omega, r, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j])
    
    diffFunc(arg, k4, params, omega, r, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j])
    
    diffFunc(arg, k5, params, omega, r, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j])
    
    diffFunc(arg, k6, params, omega, r, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j])
    
    diffFunc(arg, k7, params, omega, r, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j])
    
    diffFunc(arg, k8, params, omega, r, H)
    for j in range(dimension):
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j])


if __name__ == "__main__":
    fig, ax = plt.subplots()
    x1, x2 = [], []

    # Начальная точка из хаоса
    main_traj = np.array([-34.697048094358898, -62.831273385025732, 62.2309106274670827, 
                          -0.16059001362453112, -0.29080499936213250, -0.94320904356884894]) # M1 M2 M3 gamma1 gamma2 gamma3 
    M = np.array([-34.697048094358898, -62.831273385025732, 62.2309106274670827])
    gamma = np.array([-0.16059001362453112, -0.29080499936213250, -0.94320904356884894])

    # Параметры из хаоса
    params = [0.485, 2, 6, 7, 9, 4, 1, 752, 100] # d I1 I2 I3 a1 a2 h E g0
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    a = [a1,a2]
    h, E, g0 = params[6], params[7], params[8]
    for i in range (1):
        # print('GAMMA', gamma)
        # print("проверка геом интеграла", np.dot(gamma, gamma))

        #теперь вычислим вектор r
        r = calc_r_vec(a1, a2, h, gamma, d)
        print("R", r)

        #теперь вычислим вектор omega
        omega = calc_omega(I1, I2, I3, r, gamma, M)
        print("omega", omega)
        # m3, omega1, omega2, omega3 = calc_m3_and_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
        # M[2] = m3

        # 
        # dgamma = np.cross(gamma, omega)

        # dr = np.ndarray(3)
        # dr[0] = -0.1e1 / gamma[2] * (a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * dgamma[0] - 0.1e1 / gamma[2] * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * dgamma[1] + pow(gamma[2], -0.2e1) * ((a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1]) * dgamma[2]
        # dr[1] = -0.1e1 / gamma[2] * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * dgamma[0] - 0.1e1 / gamma[2] * (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * dgamma[1] + pow(gamma[2], -0.2e1) * ((a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1]) * dgamma[2]
        # dr[2] = (0.2e1 * (a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + 0.2e1 * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1]) * pow(gamma[2], -0.2e1) * dgamma[0] / 0.2e1 + (0.2e1 * (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + 0.2e1 * (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1]) * pow(gamma[2], -0.2e1) * dgamma[1] / 0.2e1 - (((a[0] * pow(np.cos(d), 0.2e1) + a[1] * pow(np.sin(d), 0.2e1)) * gamma[0] + (a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[1]) * gamma[0] + ((a[0] * np.cos(d) * np.sin(d) - a[1] * np.cos(d) * np.sin(d)) * gamma[0] + (a[0] * pow(np.sin(d), 0.2e1) + a[1] * pow(np.cos(d), 0.2e1)) * gamma[1]) * gamma[1]) * pow(gamma[2], -0.3e1) * dgamma[2]


        # Попробуем обратно вернуться к прежним осям
        # M = QC @ M
        # r = QC @ r
        # gamma = QC @ gamma
        # omega = QC @ omega
        # main_traj = np.array([M[0],M[1], M[2], gamma[0],gamma[1],gamma[2]])
        main_traj = np.array([*M, *gamma])
        print("проверка интеграла энергии", (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma))
        dverkStep(main_traj, 6, sys, params, 0.001, omega, r)
        M[0], M[1], M[2], gamma[0], gamma[1], gamma[2] = main_traj[0], main_traj[1], main_traj[2], main_traj[3], main_traj[4], main_traj[5]
        x1.append(i)
        # omega_y.append(omega3)
        
        omega_y.append((1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma))

    ax.plot(x1, omega_y)
    plt.show()
