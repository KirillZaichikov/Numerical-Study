import math as m
import numpy as np
import matplotlib.pyplot as plt

eps = 0.007
alpha = 0.05
beta = 0.2
p_u, p_s = 1, 1

def RevS3(ksi2, eta2, teta2):
    ####### Global
    x = 0
    for i in range(10):
        x_0 = x
        x = x_0 - (x_0+eps*m.sin(x_0)-eps*m.cos(2*x_0)-eps*m.sin(2*x_0)-teta2+ksi2+eta2+1) / \
                (1+eps*m.cos(x_0) + 2*eps*m.sin(2*x_0)-2*eps*m.cos(2*x_0))
    phi1 = x
    while (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    while (phi1 < -m.pi):
        phi1 += 2 * m.pi
    ksi1 = ksi2 - eps * m.cos(2*phi1)
    eta1 = eta2 - eps * m.sin(2*phi1)

    # print("(RevS3) phi1, ksi1, eta1=", phi1, ksi1,eta1)

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

    # while (phi0 > m.pi):
    #     phi0 = phi0 - 2 * m.pi
    # while (phi0 < -m.pi):
    #     phi0 += 2 * m.pi

    # LOCAL 2
    # ksi0 = (p_s*p_u)**((2*eps)/(alpha+eps))*ksi1*(ksi1**2+eta1**2)**(-eps/(alpha+eps))
    # eta0 = (p_s*p_u)**((2*eps)/(alpha+eps))*eta1*(ksi1**2+eta1**2)**(-eps/(alpha+eps))
    # teta0 = phi1-((beta-eps)/(alpha+eps))*np.log(p_s*p_u/(np.sqrt(ksi1**2+eta1**2)))-np.atan2(eta1,ksi1)

    while (teta0 > m.pi):
        teta0 = teta0 - 2 * m.pi
    while (teta0 < -m.pi):
        teta0 += 2 * m.pi
    
    # print("(RevS3) ksi0, eta0, teta0=", ksi0,eta0,teta0)

    return ksi0, eta0, teta0

def DirS3(ksi0, eta0, teta0):
    # print("(DirS3)  ksi0, eta0, teta0=", ksi0, eta0, teta0)
    # LOCAL 1
    # r_factor = (p_s * p_u) / m.sqrt(ksi0 * ksi0 + eta0 * eta0)
    # scale = pow(r_factor, (-2.0 * eps) / (alpha - eps))
    # Phi0 = m.atan2(eta0, ksi0)
    # angle = (2.0 * eps / (alpha - eps)) * m.log(r_factor)
        
    # ksi1 = scale * ( ksi0 * m.cos(angle) + eta0 * m.sin(angle))
    # eta1 = scale * (-ksi0 * m.sin(angle) + eta0 * m.cos(angle))
    # phi1 = teta0 + Phi0 + ((beta - eps) / (alpha - eps)) * m.log(r_factor)
    # if (phi1 > m.pi):
    #     phi1 = phi1 - 2 * m.pi
    # elif (phi1 < -m.pi):
    #     phi1 += 2 * m.pi

    # LOCAL 2

    ro = m.sqrt(ksi0**2+eta0**2)/p_s
    teta0 = teta0
    phi0 = m.atan2(eta0, ksi0) + teta0

    # print("ro, teta0, phi0")
    # print(ro, teta0, phi0)

    r = p_s * (ro / p_u) ** ((alpha+eps)/(alpha-eps))
    teta1 = ((beta + eps) / (alpha - eps)) * m.log(p_u/ro) + teta0
    phi1 = ((beta - eps) / (alpha - eps)) * m.log(p_u/ro) + phi0

    while (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    while (phi1 < -m.pi):
        phi1 += 2 * m.pi

    # print("r, teta1, phi1")
    # print(r, teta1, phi1)

    ksi1 = r * p_u * m.cos(phi1-teta1)
    eta1 = r * p_u * m.sin(phi1-teta1)
    phi1 = phi1
    
    # print("(DirS3) phi1, ksi1, eta1=", phi1, ksi1, eta1)

    # LOCAL 3
    # r_factor = (p_s * p_u)
    # scale = pow(r_factor, (-2.0 * eps) / (alpha - eps))
    # Phi0 = m.atan2(eta0, ksi0)
    # angle = (2.0 * eps / (alpha - eps)) * m.log(r_factor)
        
    # ksi1 = scale * ksi0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # eta1 = scale * eta0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # phi1 = teta0 + Phi0 + ((beta - eps) / (alpha - eps)) * m.log(r_factor/(np.sqrt(ksi0**2+eta0**2)))
    while (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    while (phi1 < -m.pi):
        phi1 += 2 * m.pi

    # GLOBAL
    ksi2 = ksi1 + eps * m.cos(2 * phi1)
    eta2 = eta1 + eps * m.sin(2 * phi1)
    teta2 = phi1 + 1 + ksi1 + eta1 + eps * m.sin(phi1)

    while (teta2 > m.pi):
        teta2 -= 2 * m.pi
    while (teta2 < -m.pi):
        teta2 += 2 * m.pi

    return ksi2, eta2, teta2


iter_num=1
y1, y2 = [], []
ksi0, eta0, teta0 =  -0.02606655, -0.04264003, -0.10349153
# old = array([-0.02606655, -0.04264003, -0.10349153])
# ksi0, eta0, teta0 = DirS3(ksi0,eta0,teta0)
# print(ksi0,eta0,teta0)

# ksi0, eta0, teta0 = RevS3(ksi0,eta0,teta0)
# print(ksi0,eta0,teta0)

for i in range(iter_num):
    y1.append([ksi0, eta0, teta0])
    ksi0, eta0, teta0 = RevS3(ksi0,eta0,teta0)
    print([ksi0, eta0, teta0])
y1.append([ksi0, eta0, teta0])
print(len(y1))
# print("rev start")
for i in range(iter_num):
    y2.append([ksi0, eta0, teta0])
    ksi0, eta0, teta0 = DirS3(ksi0,eta0,teta0)
    print([ksi0, eta0, teta0])
# y2.append([ksi0, eta0, teta0])
# print(len(y2))
# print(y1, "\n",y2)
# for i in range(iter_num):
#     print("rev", y2[i])
#     print("dir", y1[-1 - i])

t = np.arange(0, iter_num+1)
yn1 = np.array(y1)
yn2 = np.array(y2)
y1 = yn1.T
y2 = yn2.T


fig, ax = plt.subplots()
ax.plot(t, y1[2])
# ax.plot(t, y2[2][::-1])
fig1, ax1 = plt.subplots()
ax1.plot(y1[1],y1[2],linestyle="", marker = '.', ms=5)
# ax1.plot(y2[1],y2[2],linestyle="", marker = '.', ms=5)
plt.show()