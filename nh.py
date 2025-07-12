import math as m
import numpy as np
from utils.integrator import *
import matplotlib.pyplot as plt

omega_y = []
# def sys(res, m1, m2, m3, gamma1, gamma2, gamma3, d, I1, I2, I3, a1, a2, h, E, g0): # масса 1, gamma3 считаем из геом интеграла


def calc_gamma3(gamma1, gamma2):
    # из геом интеграла g1^2 + g2^2 + g3^2 = 1 выразим gamma3 пусть с плюсом
    print("проверка геометрического интеграла", gamma1**2+gamma2**2+m.sqrt(1-gamma1**2-gamma2**2)**2)
    return m.sqrt(1-gamma1**2-gamma2**2)

def calc_r_vec(a1, a2, h, gamma):
    r1 = -(a1 * gamma[0]) / gamma[2]
    r2 = -(a2 * gamma[1]) / gamma[2]
    r3 = -h + (1/2)*(a1 * gamma[1] ** 2 + a2 * gamma[1] ** 2)/(gamma[2] ** 2)
    return r1, r2, r3

def calc_m3_and_omega(d, I1, I2, I3, a1, a2,h,E,g0, r, gamma):
    r1,r2,r3 = r[0], r[1], r[2]
    # gamma1, gamma2, gamma3 = gamma[0], gamma[1], gamma[2]
    A = np.array([[I1 * (m.cos(d))**2 + I2*(m.sin(d))**2 + r2 ** 2 + r3 ** 2, 
                  (I1-I2) * m.cos(d) * m.sin(d) - r1 * r2,
                  - r1 * r2],
                 [(I1-I2) * m.cos(d) * m.sin(d) - r1 * r2,
                  I1 * (m.sin(d))**2 + I2*(m.cos(d))**2 + r1 ** 2 + r3 ** 2,
                  - r2 * r3],
                 [- r1 * r3, - r2 * r3, I3 + r1 ** 2 + r2 ** 2]])
    A_inv = np.linalg.inv(A)
    # print('проверка обратной матрицы', A @ A_inv)
    a11, a12, a13 = A_inv[0][0], A_inv[0][1], A_inv[0][2]
    a21, a22, a23 = A_inv[1][0], A_inv[1][1], A_inv[1][2]
    a31, a32, a33 = A_inv[2][0], A_inv[2][1], A_inv[2][2]

    B = a13 * m1 + a31 * m1 + a23 * m2 + a32 * m2
    C = a11 * m1 ** 2 + (a12+a21)*m1*m2 + a22 * m2 ** 2
    D = (B/2)**2-4*a33/2*(C/2-g0 * np.dot(r, gamma) - E)
    m3 = (-B/2 + m.sqrt(D)) / a33                                              # взял m3 тот что с плюсом
    M = np.array([m1,m2,m3])
    omega = A_inv @ np.array([m1, m2, m3])
    omega1 = omega[0]
    omega2 = omega[1]
    omega3 = omega[2]
    omega = np.array([omega1,omega2,omega3])
    # omega_y.append(omega1)
    #проверка
    print("проверка интеграла энергии", 1/2 * np.dot(M, omega) - g0 * np.dot(r, gamma)-E)
    return m3, omega1, omega2, omega3

def sys(state, res, params, omega, H): # масса 1, gamma3 считаем из геом интеграла
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]
    m1, m2, gamma1, gamma2 = state[0], state[1], state[2], state[3]
    omega1, omega2, omega3 = omega[0], omega[1], omega[2]
    # # из геом интеграла g1^2 + g2^2 + g3^2 = 1 выразим gamma3 пусть с плюсом
    # gamma3 = m.sqrt(1-gamma1**2-gamma2**2)

    # соотношение между r и gamma, 
    # a1, a2, h -- геометрические параметры 
    # r1 = -(a1 * gamma1) / gamma3
    # r2 = -(a2 * gamma2) / gamma3
    # r3 = -h + (1/2)*(a1 * gamma1 ** 2 + a2 * gamma2 ** 2)/(gamma3 ** 2)
    # r = np.array([r1,r2,r3])

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

    # Вычисление p и omegar 
    r1_ = - a1 * ((gamma2 * omega3 - gamma3 * omega2) * gamma3 - 
                  (gamma1 * omega2 - gamma2*omega1)*gamma1)/ (gamma3**2) 
    r2_ = - a2 * ((gamma3 * omega1 - gamma1 * omega3) * gamma3 - 
                  (gamma1 * omega2 - gamma2*omega1)*gamma2)/ (gamma3**2) 
    r3_ = (a1 * gamma1 * (gamma2 * omega3 - gamma3 * omega2) + 
           a2 * gamma2 * (gamma3 * omega1 - gamma1 * omega3)) / (gamma3**2) - \
          (a1 * gamma1 ** 2 + a2 * gamma2 ** 2)*(gamma1 * omega2 - gamma2 * omega1)/(gamma3 ** 3)
    p = r1 * r1_ + r2 * r2_ + r2 * r3_
    omegar = r1_ * omega1 + r2_ * omega2 + r3_ * omega3

    m1_ = m2 * omega3 - m3 * omega2 + omega1 * p - r1 * omegar + g0 * (r2 * gamma3 - gamma2 * r3)
    m2_ = m3 * omega1 - m1 * omega3 + omega2 * p - r2 * omegar + g0 * (r3 * gamma1 - gamma3 * r1)
    # m3_ = m1 * omega2 - m2 * omega1 + omega3 * p - r3 * omegar + g0 * (r1 * gamma2 - gamma1 * r2)
    gamma1_ = gamma2 * omega3 - gamma3 * omega2
    gamma2_ = gamma3 * omega1 - gamma1 * omega3
    # gamma3_ = gamma1 * omega2 - gamma2 * omega1

    res[0], res[1], res[2], res[3]= m1_, m2_, gamma1_, gamma2_


def dverkStep(val, dimension, diffFunc, params, step, omega, H=1) -> None:
    k1, k2, k3, k4, k5, k6, k7, k8 = np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), \
        np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension)
    arg = np.zeros(dimension)
    diffFunc(val, k1, params, omega, H)
    for j in range(dimension):
        arg[j] = val[j] + step / 6. * k1[j]
    
    diffFunc(arg, k2, params, omega,H)
    for j in range(dimension):
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j])
    
    diffFunc(arg, k3, params, omega,H)
    for j in range(dimension):
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j])
    
    diffFunc(arg, k4, params, omega,H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j])
    
    diffFunc(arg, k5, params, omega,H)
    for j in range(dimension):
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j])
    
    diffFunc(arg, k6, params, omega,H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j])
    
    diffFunc(arg, k7, params, omega,H)
    for j in range(dimension):
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j])
    
    diffFunc(arg, k8, params, omega,H)
    for j in range(dimension):
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j])


if __name__ == "__main__":
    fig, ax = plt.subplots()
    x1, x2 = [], []

    main_traj = np.array([1, 1, 0.5, 0.5]) # M1 M2 M3 gamma1 gamma2 gamma3 Сейчас gamma3 можно ставить любое

    params = [0.2, 5, 6, 7, 9, 4, 1, 1380, 100] # d I1 I2 I3 a1 a2 h E g0
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    for i in range (100):
        m1, m2, gamma1, gamma2 = main_traj[0], main_traj[1], main_traj[2], main_traj[3]

        gamma3 = calc_gamma3(gamma1, gamma2)
        gamma = np.array([gamma1,gamma2,gamma3])

        r1, r2, r3 = calc_r_vec(a1, a2, h, gamma)
        r = np.array([r1,r2,r3])

        m3, omega1, omega2, omega3 = calc_m3_and_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma)
        omega = np.array([omega1, omega2, omega3])
        M = np.array([m1, m2, m3])

        dverkStep(main_traj, 4, sys, params, 0.01, omega)
        x1.append(i)
        omega_y.append(omega3)
    ax.plot(x1, omega_y)
    plt.show()
