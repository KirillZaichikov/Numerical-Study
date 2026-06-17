import numpy as np
import math as m
# from numba import jit

def lerFRM_3D_map(state, res, params):
    eps, alpha, beta, p = params
    # print(eps, alpha, p_s, p_u)
    ksi0, eta0, teta0 = state

    # LOCAL 1
    # r_factor = (p**2) / np.sqrt(ksi0**2 + eta0**2)
    # scale = (r_factor) ** ((-2 * eps) / (alpha - eps))

    # Phi0 = np.arctan2(eta0, ksi0)
    # # Phi1 = (-(2 * eps) / (alpha - eps)) * np.log(r_factor) + Phi0
    # # if Phi1>np.pi:
    # #     Phi1 -= 2 * np.pi
    # # elif Phi1<-np.pi:
    # #     Phi1 += 2 * np.pi

    # # ksi1 = p_s * p_u * scale * np.cos(Phi1)
    # # eta1 = p_s * p_u * scale * np.sin(Phi1)
    # # phi1 = Phi1 - teta0
    # ksi1 = scale * ( ksi0 * np.cos( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.sin(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # eta1 = scale * (-ksi0 * np.sin( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.cos(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # # phi1 = (teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)) % (2*np.pi)
    # phi1 =  teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)

    # LOCAL 2
    ro = np.sqrt(ksi0**2+eta0**2)/p
    teta0 = teta0
    phi0 = np.atan2(eta0, ksi0) + teta0

    # while (phi0 > m.pi):
    #     phi0 -= 2 * m.pi
    # while (phi0 < -m.pi):
    #     phi0 += 2 * m.pi

    # print("ro, teta0, phi0")
    # print(ro, teta0, phi0)

    r = p * (ro / p) ** ((alpha+eps)/(alpha-eps))
    teta1 = ((beta + eps) / (alpha - eps)) * np.log(p/ro) + teta0
    phi1 = (((beta - eps) / (alpha - eps)) * np.log(p/ro) + phi0)

    # print(r, teta1, phi1)

    while (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    while (phi1 < -m.pi):
        phi1 += 2 * m.pi

    # while (teta1 > m.pi):
    #     teta1 = teta1 - 2 * m.pi
    # while (teta1 < -m.pi):
    #     teta1 += 2 * m.pi
    # print("r, teta1, phi1")
    # print(r, teta1, phi1)
    ksi1 = r * p * np.cos(phi1-teta1)
    eta1 = r * p * np.sin(phi1-teta1)
    phi1 = phi1

    # print(ksi1, eta1, phi1)
    # LOCAL 3
    # r_factor = (p * p)
    # scale = pow(r_factor, (-2.0 * eps) / (alpha - eps))
    # Phi0 = m.atan2(eta0, ksi0)
    # # angle = (2.0 * eps / (alpha - eps)) * m.log(r_factor)
        
    # ksi1 = scale * ksi0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # eta1 = scale * eta0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # phi1 = teta0 + Phi0 + ((beta - eps) / (alpha - eps)) * m.log(r_factor/(np.sqrt(ksi0**2+eta0**2)))
    # while (phi1 > m.pi):
    #     phi1 = phi1 - 2 * m.pi
    # while (phi1 < -m.pi):
    #     phi1 += 2 * m.pi

    # ******** S3
    # ksi2 = ksi1 + eps * np.cos(2*phi1)
    # eta2 = eta1 + eps * np.sin(2*phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    # teta2 = phi1

    # ******** S2
    ksi2 = ksi1 + eps * np.cos(phi1)
    eta2 = eta1 + eps * np.sin(phi1)
    # # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S1
    # ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    # eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    # teta2 =  phi1

    # print(ksi2, eta2, teta2)
    while teta2>np.pi:
        teta2 -= 2 * np.pi
    while teta2<-np.pi:
        teta2 += 2 * np.pi

    res[0] = ksi2 
    res[1] = eta2
    res[2] = teta2