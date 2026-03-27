import math as m
import numpy as np
import matplotlib.pyplot as plt

eps = 0.1
alpha = 0.8
beta = 0.2
p_u, p_s = 1, 1

def RevS3(ksi2, eta2, teta2):
    ####### Global S1
    x = 2.5
    for i in range(10):
        x_0 = x
        x = x_0 - (x_0+0.25*eps*m.sin(x_0)-0.25*eps*m.cos(x_0)-teta2+ksi2+eta2-0.5*eps+1) / \
                (1+0.25*eps*m.cos(x_0) + 0.25*eps*m.sin(x_0))
    phi1 = x
    if (phi1 > m.pi):
        phi1 -= 2 * m.pi
    elif (phi1 < -m.pi):
        phi1 += 2 * m.pi
    ksi1 = ksi2 - eps * (0.5 + 0.25 * m.cos(phi1))
    eta1 = eta2 - 0.75 * eps * m.sin(phi1)

    ####### Local 1
    eps_ = -eps
    alpha_ = -alpha
    beta_ = -beta

    r = m.sqrt(ksi1**2 + eta1**2) / p_u
    teta1 = phi1 - m.atan2(eta1, ksi1)
    while (teta1 > m.pi):
        teta1 = teta1 - 2 * m.pi
    while (teta1 < -m.pi):
        teta1 += 2 * m.pi
    phi1 = phi1

    ro = p_u * (r / p_s) ** ((alpha_-eps_)/(alpha_+eps_))
    teta0 = ((beta_ + eps_) / (-(alpha_ + eps_))) * m.log(p_s/r) + teta1
    phi0 =  ((beta_ - eps_) / (-(alpha_ + eps_))) * m.log(p_s/r) + phi1

    ksi0 = ro * p_s * m.cos(phi0-teta0)
    eta0 = ro * p_s * m.sin(phi0-teta0)
    phi0 = phi0

    while (teta0 > m.pi):
        teta0 = teta0 - 2 * m.pi
    while (teta0 < -m.pi):
        teta0 += 2 * m.pi

    return ksi0, eta0, teta0


def DirS3(ksi0, eta0, teta0):
    print("(DirS3)  ksi0, eta0, teta0=", ksi0, eta0, teta0)
    #LOCAL 1
    # r_factor = (p_s * p_u) / m.sqrt(ksi0 * ksi0 + eta0 * eta0)
    # scale = pow(r_factor, (-2.0 * alpha) / (alpha - eps))
    # Phi0 = m.atan2(eta0, ksi0)
    # angle = (2.0 * eps / (alpha - eps)) * m.log(r_factor)
        
    # ksi1 = scale * (ksi0 * m.cos(angle) + eta0 * m.sin(angle))
    # eta1 = scale * (-ksi0 * m.sin(angle) + eta0 * m.cos(angle))
    # phi1 = teta0 + Phi0 + ((beta - eps) / (alpha - eps)) * m.log(r_factor)
    # if (phi1 > m.pi):
    #     phi1 = phi1 - 2 * m.pi
    # elif (phi1 < -m.pi):
    #     phi1 += 2 * m.pi

    #LOCAL 2

    ro = m.sqrt(ksi0**2+eta0**2)/p_s
    teta0 = teta0
    phi0 = m.atan2(eta0, ksi0) + teta0

    print("ro, teta0, phi0")
    print(ro, teta0, phi0)

    r = p_s * (ro / p_u) ** ((alpha+eps)/(alpha-eps))
    teta1 = ((beta + eps) / (alpha - eps)) * m.log(p_u/ro) + teta0
    phi1 = ((beta - eps) / (alpha - eps)) * m.log(p_u/ro) + phi0

    if (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    elif (phi1 < -m.pi):
        phi1 += 2 * m.pi

    print("r, teta1, phi1")
    print(r, teta1, phi1)

    ksi1 = r * p_u * m.cos(phi1-teta1)
    eta1 = r * p_u * m.sin(phi1-teta1)
    phi1 = phi1
    
    print("(DirS3) phi1, ksi1, eta1=", phi1, ksi1, eta1)
    # GLOBAL
    ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    if (teta2 > m.pi):
        teta2 -= 2 * m.pi
    elif (teta2 < -m.pi):
        teta2 += 2 * m.pi

    return ksi2, eta2, teta2


iter_num=100
y1, y2 = [], []
ksi0, eta0, teta0 = 0.1, 0.001, 0.2
for i in range(iter_num):
    ksi0, eta0, teta0 = DirS3(ksi0,eta0,teta0)
    print([ksi0, eta0, teta0])
    y1.append([ksi0, eta0, teta0])
print("rev start")
for i in range(iter_num):
    y2.append([ksi0, eta0, teta0])
    print([ksi0, eta0, teta0])
    ksi0, eta0, teta0 = RevS3(ksi0,eta0,teta0)

t = np.arange(0, iter_num)
yn1 = np.array(y1)
yn2 = np.array(y2)
y1 = yn1.T
y2 = yn2.T

fig, ax = plt.subplots()
ax.plot(t, y1[2])
ax.plot(t, y2[2][::-1])
fig1, ax1 = plt.subplots()
ax1.plot(y1[1],y1[2],linestyle="", marker = '.', ms=5)
ax1.plot(y2[1],y2[2],linestyle="", marker = '.', ms=5)
plt.show()
