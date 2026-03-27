import math as m
import numpy as np

eps = 0.05
alpha = 0.8
beta = 0.2

p_s = 1
p_u = 1

ksi0 = 0.01
eta0 = 0.1
teta0 = 0.001

print("ksi0, eta0, teta0")
print(ksi0, eta0, teta0)

# ro = m.sqrt(ksi0**2+eta0**2)/p_s
# teta0 = teta0
# phi0 = m.atan2(eta0, ksi0) + teta0

# print("ro, teta0, phi0")
# print(ro, teta0, phi0)

# r = p_s * (ro / p_u) ** ((alpha+eps)/(alpha-eps))
# teta1 = ((beta + eps) / (alpha - eps)) * m.log(p_u/ro) + teta0
# phi1 = ((beta - eps) / (alpha - eps)) * m.log(p_u/ro) + phi0

# print("r, teta1, phi1")
# print(r, teta1, phi1)

# ksi1 = r * p_u * m.cos(phi1-teta1)
# eta1 = r * p_u * m.sin(phi1-teta1)
# phi1 = phi1

# LOCAL 2
r_factor = (p_s * p_u) / m.sqrt(ksi0 * ksi0 + eta0 * eta0)
scale = pow(r_factor, (-2.0 * alpha) / (alpha - eps))
Phi0 = m.atan2(eta0, ksi0)
angle = (2.0 * eps / (alpha - eps)) * m.log(r_factor)
    
ksi1 = scale * (ksi0 * m.cos(angle) + eta0 * m.sin(angle))
eta1 = scale * (-ksi0 * m.sin(angle) + eta0 * m.cos(angle))
phi1 = teta0 + Phi0 + ((beta - eps) / (alpha - eps)) * m.log(r_factor)

print("ksi1, eta1, phi1")
print(ksi1, eta1, phi1)

###### В ОБРАТНОМ ВРЕМЕНИ
# eps = -0.05
# alpha = -0.8
# beta = -0.2

# r = m.sqrt(ksi1**2 + eta1**2) / p_u
# teta1 = phi1 - m.atan2(eta1, ksi1)
# phi1 = phi1

# print("r, teta1, phi1")
# print(r, teta1, phi1)

# ro = p_u * (r / p_s) ** ((alpha-eps)/(alpha+eps))
# teta0 = ((beta + eps) / (-(alpha + eps))) * m.log(p_s/r) + teta1
# phi0 =  ((beta - eps) / (-(alpha + eps))) * m.log(p_s/r) + phi1

# print("ro, teta0, phi0")
# print(ro, teta0, phi0)

# ksi0 = ro * p_s * m.cos(phi0-teta0)
# eta0 = ro * p_s * m.sin(phi0-teta0)
# phi0 = phi0

# LOCAL 2
ksi0 = (p_s*p_u)**((2*eps)/(alpha+eps))*ksi1*(ksi1**2+eta1**2)**(-eps/(alpha+eps))
eta0 = (p_s*p_u)**((2*eps)/(alpha+eps))*eta1*(ksi1**2+eta1**2)**(-eps/(alpha+eps))
teta0 = phi1-((beta-eps)/(alpha+eps))*np.log((p_s*p_u)/(np.sqrt(ksi1**2+eta1**2)))+np.atan2(eta1,ksi1)


print("ksi0, eta0, teta0")
print(ksi0, eta0, teta0)
