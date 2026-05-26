'''
тут лежат функции для координат Кузнецова
'''

import numpy as np

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

def sysRev(state, res, params, H): # масса 1
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
    
    dgamma = -np.cross(gamma, omega)
    dgamma_tmp = np.cross(gamma, omega)
    Jr = np.array([[-a1/gamma[2], 0, a1 * gamma[0]/(gamma[2]**2)],
                   [0, -a2/gamma[2], a2 * gamma[1]/(gamma[2]**2)],
                   [a1*gamma[0]/(gamma[2]**2), a2*gamma[1]/(gamma[2]**2), -(a1*gamma[0]**2 + a2 * gamma[1]**2)/(gamma[2]**3)]])
    dr = Jr @ dgamma_tmp

    dM = -(np.cross(M, omega) + np.cross(dr, np.cross(omega, r)) + g0 * np.cross(r, gamma))

    res[0], res[1], res[2], res[3], res[4], res[5] = dM[0], dM[1], dM[2], dgamma[0], dgamma[1], dgamma[2]
    # res[0], res[1], res[2], res[3], res[4], res[5] = m1_, m2_, m3_, gamma1_, gamma2_, gamma3_

def calc_integrals(M, gamma, d, I1, I2, I3, a1, a2, h, E, g0):
    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    En = (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma)
    mg = np.dot(gamma, gamma)
    print(f"Energy is {En}", f'geom int {mg}', sep="\n")
    return En, mg
