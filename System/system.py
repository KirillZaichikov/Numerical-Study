import math as m
import numpy as np
from numba import jit

@jit(nopython=True, cache=True)
def ctg(x): return 1 / m.tan(x)

@jit(nopython=True, cache=True)
def main_sys(start_point, res, params, H):
    ro, mpar, i1, i2, a0, g, en = params
    ksi, phi, teta = start_point

    # print(params)
    # print(start_point)
    omega2 = 2*(en - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
            (i1 * np.cos(ksi-phi) ** 2 + i2 * np.sin(ksi-phi) ** 2 + \
             mpar * np.cos(teta)**2*(ro*np.cos(ksi)-a0 * np.sin(ksi-phi))**2)
    # print(omega2)
    omega = np.sqrt(omega2)

    if omega2<0: print("OMEGA less than 0")

    N = ro * omega2 * (m.cos(ksi)**2-m.cos(teta)**2) - \
        m.sin(teta)*(g - a0 * omega2 * m.sin(teta) * m.sin(ksi-phi) * m.cos(ksi))

    D = i1*i2 + mpar * m.cos(teta)**2 * (i1 * (a0 + ro * m.sin(phi)) ** 2 + \
                                      i2 * ro ** 2 * m.cos(phi) ** 2)
    
    # Сама система
    ksi1 = -ctg(teta) * (omega * m.sin(ksi) + \
                         ((N*mpar)/(D*omega)) * (i1*(a0+ro*m.sin(phi))*m.cos(ksi-phi) + \
                                                 i2*ro*m.cos(phi)*m.sin(ksi-phi)) )
    phi1 = -omega * ctg(teta) * m.sin(ksi)
    teta1 = omega * m.cos(ksi)

    res[0] = ksi1
    res[1] = phi1
    res[2] = teta1

@jit(nopython=True, cache=True)
def main_sys_poincare(start_point, res, params, H):
    ro, mpar, i1, i2, a0, g, eps = params
    ksi, phi, teta = start_point

    # print(params)
    # print(start_point)
    omega2 = 2*(eps - mpar * g * m.sin(teta) * (ro + a0*m.sin(phi)) ) / \
            (i1 * m.cos(ksi-phi) ** 2 + i2 * m.sin(ksi-phi) ** 2 + \
             mpar * m.cos(teta)**2*(ro*m.cos(ksi)-a0 * m.sin(ksi-phi))**2)
    # print(omega2)
    omega = m.sqrt(omega2)
    if omega2<0: print("OMEGA less than 0")

    N = ro * omega2 * (m.cos(ksi)**2-m.cos(teta)**2) - \
        m.sin(teta)*(g - a0 * omega2 * m.sin(teta) * m.sin(ksi-phi) * m.cos(ksi))

    D = i1*i2 + mpar * m.cos(teta)**2 * (i1 * (a0 + ro * m.sin(phi)) ** 2 + \
                                      i2 * ro ** 2 * m.cos(phi) ** 2)
    
    # Сама система
    ksi1 = -ctg(teta) * (omega * m.sin(ksi) + \
                         ((N*mpar)/(D*omega)) * (i1*(a0+ro*m.sin(phi))*m.cos(ksi-phi) + \
                                              i2*ro*m.cos(phi)*m.sin(ksi-phi)) )
    phi1 = -omega * ctg(teta) * m.sin(ksi)
    teta1 = omega * m.cos(ksi)

    res[0] = ksi1  / (- omega * ctg(teta) * m.sin(ksi))
    res[1] = phi1  / (- omega * ctg(teta) * m.sin(ksi))
    res[2] = teta1 / (- omega * ctg(teta) * m.sin(ksi))
    # res[0] = ksi1 / H
    # res[1] = phi1 / H
    # res[2] = teta1 / H
    
@jit(nopython=True, cache=True)
def omega_sys(start_point, res, params, H):
    ro, mpar, i1, i2, a0, g, en = params
    ksi, phi, teta, omega = start_point

    N = ro * omega**2 * (np.cos(phi+ksi)**2 - np.cos(teta)**2) - \
        np.sin(teta) * (g - a0 * omega ** 2 * np.sin(teta) * np.sin(ksi)*np.cos(phi+ksi))

    # Сама система
    ksi1 = -ctg(teta) * (omega * m.sin(ksi) + \
                         ((N*mpar)/(D*omega)) * (i1*(a0+ro*m.sin(phi))*m.cos(ksi-phi) + \
                                                 i2*ro*m.cos(phi)*m.sin(ksi-phi)) )
    phi1 = -omega * ctg(teta) * m.sin(ksi)
    teta1 = omega * m.cos(ksi)

    res[0] = ksi1
    res[1] = phi1
    res[2] = teta1

@jit(nopython=True, cache=True)
def omega12_sys(start_point, res, params, H):
    ro, mpar, i1, i2, a0, g, en = params
    dzeta, phi, teta = start_point

    omega2 = 2*(en - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
        (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
            mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2)
    if omega2<0: print("OMEGA less than 0")
    omega = np.sqrt(omega2)

    omega_1 = omega * np.cos(dzeta)
    omega_2 = omega * np.sin(dzeta)


    
    # Сама система
    ksi1 = -ctg(teta) * (omega * m.sin(ksi) + \
                         ((N*mpar)/(D*omega)) * (i1*(a0+ro*m.sin(phi))*m.cos(ksi-phi) + \
                                                 i2*ro*m.cos(phi)*m.sin(ksi-phi)) )
    phi1 = -omega * ctg(teta) * m.sin(ksi)
    teta1 = omega * m.cos(ksi)

    res[0] = ksi1
    res[1] = phi1
    res[2] = teta1

@jit(nopython=True, cache=True)
def omegaGamma_sys(start_point, res, params, H):
    omega1, omega2, gamma1, gamma2, gamma3 = start_point
    rho, mpar, I_1, I_2, a0, g, en = params

    w_1 = omega1
    w_2 = omega2
    _g1 = gamma1
    _g2 = gamma2
    _g3 = gamma3

    __phi = np.atan2(_g1, _g2)
    
    coun=0
    while ((__phi<0.0) or (__phi>2.0*np.pi)):
        if (__phi<0.0):
            __phi += 2.0 * np.pi
        else:
            __phi -= 2.0 * np.pi
        coun+=1
        if coun>50:
            print("ERRORERRORERROR")
            break
    
    if __phi>=2*np.pi or __phi<=0: print("ERROR (mod phi):",__phi)

    t1 = _g3 * _g3
    t2 = 1 - t1
    if t2 <0: print("ERROR (t2 less than 0)", t2)
    t3 = np.sqrt(t2)
    t4 = -_g1*w_2 + _g2*w_1
    t41 = -np.sin(__phi)*w_2 + np.cos(__phi)*w_1
    t5 = w_1*w_1
    t6 = w_2*w_2
    t7 = t5 + t6
    n1 = -(w_2*t4*a0 + g)*t3
    n2 = -rho*t7*t1
    n3 = t41*t41*rho
    br = a0*t3 + _g1*rho
    t21 = rho*rho
    t8 = mpar*rho*_g2*I_2
    t9 = mpar*br*I_1
    t10 = _g1*_g1
    t11 = _g2*_g2
    t12 = a0*a0
    t13 = t1*t1
    dnm = 2 * I_1*t3*_g1*t1*a0*mpar*rho - I_1*t13*t12*mpar + (mpar*(I_1*t10 + I_2*t11)*t21 + I_1*(t12*mpar - I_2))*t1 + I_1*I_2

    res[0] = t8*_g3*(n1 + n2 + n3) / dnm
    res[1] = -t9*_g3*(n1 + n2 + n3) / dnm

    res[2] = -w_2*_g3
    res[3] = w_1*_g3
    res[4] = (_g1*w_2 - _g2*w_1)

@jit(nopython=True, cache=True)
def omegaGamma_sys_poincare(start_point, res, params, H):
    omega1, omega2, gamma1, gamma2, gamma3 = start_point
    rho, mpar, I_1, I_2, a0, g, en = params

    w_1 = omega1
    w_2 = omega2
    _g1 = gamma1
    _g2 = gamma2
    _g3 = gamma3

    __phi = np.atan2(_g1, _g2)
    
    coun=0
    while ((__phi<0.0) or (__phi>2.0*np.pi)):
        if (__phi<0.0):
            __phi += 2.0 * np.pi
        else:
            __phi -= 2.0 * np.pi
        coun+=1
        if coun>50:
            print("ERRORERRORERROR")
            break
    
    if __phi>=2*np.pi or __phi<=0: print("ERROR (mod phi):",__phi)

    t1 = _g3 * _g3
    t2 = 1 - t1
    if t2 <0: print("ERROR (t2 less than 0)", t2)
    t3 = np.sqrt(t2)
    t4 = -_g1*w_2 + _g2*w_1
    t41 = -np.sin(__phi)*w_2 + np.cos(__phi)*w_1
    t5 = w_1*w_1
    t6 = w_2*w_2
    t7 = t5 + t6
    n1 = -(w_2*t4*a0 + g)*t3
    n2 = -rho*t7*t1
    n3 = t41*t41*rho
    br = a0*t3 + _g1*rho
    t21 = rho*rho
    t8 = mpar*rho*_g2*I_2
    t9 = mpar*br*I_1
    t10 = _g1*_g1
    t11 = _g2*_g2
    t12 = a0*a0
    t13 = t1*t1
    dnm = 2 * I_1*t3*_g1*t1*a0*mpar*rho - I_1*t13*t12*mpar + (mpar*(I_1*t10 + I_2*t11)*t21 + I_1*(t12*mpar - I_2))*t1 + I_1*I_2

    res[0] = (t8*_g3*(n1 + n2 + n3) / dnm ) / (t8*_g3*(n1 + n2 + n3) / dnm)
    res[1] = (-t9*_g3*(n1 + n2 + n3) / dnm) / (t8*_g3*(n1 + n2 + n3) / dnm)

    res[2] = -w_2*_g3 / (t8*_g3*(n1 + n2 + n3) / dnm)
    res[3] = w_1*_g3 / (t8*_g3*(n1 + n2 + n3) / dnm)
    res[4] = (_g1*w_2 - _g2*w_1) / (t8*_g3*(n1 + n2 + n3) / dnm)