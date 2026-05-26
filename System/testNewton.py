import numpy as np
import matplotlib.pyplot as plt
from params import *

def get_omega(start_point, params):
    ro, mpar, i1, i2, a0, g, eps = params
    dzeta, phi, teta = start_point
    omega2 = 2*(eps - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
        (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
            mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2)
    # print(omega2)
    omega = m.sqrt(omega2)
    return omega

def visualise(start_point, params):
    '''
    omega1,omega2,phi -> OmegaGamma
    '''
    omega1, omega2, phi = start_point
    ro, mpar, i1, i2, a0, g, en = params

    ksi = np.atan2(omega2, omega1)
    dzeta = ksi + phi
    omega = np.sqrt(omega1 ** 2 + omega2 ** 2)

    # Результат зависит от начального приближения
    teta0_ = np.linspace(0, np.pi, 1000)
    teta_x = []
    teta_y = []
    for teta0 in teta0_:
        fx0 = ((omega ** 2) / 2) * (i1*np.cos(ksi)**2+i2*np.sin(ksi)**2) + \
            (mpar/2) * np.cos(teta0) ** 2 * (omega**2) * (ro * np.cos(phi+ksi) - \
                                                          a0*np.sin(ksi))**2 + \
            mpar*g*np.sin(teta0) * (ro + a0*np.sin(phi)) - en
        dteta = teta0 + 1e-7
        fdx0 = ((((omega ** 2) / 2) * (i1*np.cos(ksi)**2+i2*np.sin(ksi)**2) + \
            (mpar/2) * np.cos(dteta) ** 2 *(omega**2) * (ro * np.cos(phi+ksi) - \
                                                        a0*np.sin(ksi))**2 + \
            mpar*g*np.sin(dteta) * (ro + a0*np.sin(phi)) - en) - fx0) / 1e-7
        teta1 = teta0 - fx0/fdx0
        teta_x.append(teta0)
        teta_y.append(teta1)

    fig,ax = plt.subplots()
    ax.set_xlim(0,np.pi)
    ax.set_ylim(0,np.pi)
    ax.plot(teta_x, teta_y)
    ax.plot(teta0_,teta0_)

    # Результат зависит от начального приближения
    teta0 = np.pi-0.1
    # teta0 = 0
    teta_traj_x = []
    teta_traj_y = []
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
        teta_traj_x.append(teta0)
        teta_traj_y.append(teta1)
        teta0 = teta1

    print("teta_traj_x", teta_traj_x)
    print("result teta", teta0)
    ax.scatter(teta_traj_x, teta_traj_y)

    gamma1=np.sin(teta0)*np.sin(phi)
    gamma2=np.sin(teta0)*np.cos(phi)
    gamma3=np.cos(teta0)

    print("result_gamma", gamma1, gamma2, gamma3)

    plt.show()

initial_point = np.array([0.01,0.5,0.01])
dzeta, phi, teta = initial_point
omega = get_omega(initial_point, params)
new_initial = np.array([omega * np.cos(dzeta-phi), omega * np.sin(dzeta-phi), phi])
visualise(new_initial, params)

