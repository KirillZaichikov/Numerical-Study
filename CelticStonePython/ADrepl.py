'''
замена из L/G H/G l в M gamma
'''


import math as m
from kuz import calc_r_vec, calc_omega, calcIntegrals
import numpy as np

def calcIntegrals(M, g0, gamma):
    r = calc_r_vec(a1, a2, h, gamma, d)
    omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
    print("(calcIntegrals) r", r, "omega", omega)
    return (1/2) * np.dot(M, omega) - g0 * np.dot(r, gamma), np.dot(gamma, gamma)

#         d     I1 I2 I3 a1 a2  h   E    g0
params = [0.485, 2, 6, 7, 9, 4, 1, 752, 100] 
d = params[0]
I1, I2, I3 = params[1], params[2], params[3]
a1, a2 = params[4], params[5]
h, E, g0 = params[6], params[7], params[8]
# первая норм точка 3.6225358 0.65189173 -0.35839377
# l = 3.6684695
# LG =  0.73725074
# HG = -0.47236776
# 3.6507227868, 0.6691444994, -0.3847016971  - Это неподвижная точка 0.485 752
l = 3.721
LG =  0.6691441646
HG = -0.3847013110

m1 = (1-LG**2)**(1/2) * m.sin(l)
m2 = (1-LG**2)**(1/2) * m.cos(l)
m3 = LG
M = np.array([m1,m2,m3]) # При такой замене вектор M будет единичной длины

gamma1 = (HG * (1-LG**2)**(1/2) + LG*(1-HG**2)**(1/2)) * m.sin(l)
gamma2 = (HG * (1-LG**2)**(1/2) + LG*(1-HG**2)**(1/2)) * m.cos(l)
gamma3 = HG * LG - (1-LG**2)**(1/2) * (1-HG**2)**(1/2) * 1
print("точка до преобразорваний", m1, m2, m3, gamma1, gamma2, gamma3)
gamma = np.array([gamma1,gamma2,gamma3]) # гамма тоже единичный
gamma = -gamma # надо повернуть
print(M, gamma)
# так получаются координаты как в хаосе, нужны как у кузнецова
QC = np.array([[np.cos(d), np.sin(d), 0],
                [-np.sin(d),np.cos(d),0],
                [0,0,1]])
M = QC @ M
gamma = QC @ gamma

# Приведем к необходимому уровню энергии и на этом все
r = calc_r_vec(a1, a2, h, gamma, d)
omega = calc_omega(d, I1, I2, I3, a1, a2, h, E, g0, r, gamma, M)
print("rg", r, omega)
tmp = M.copy()
for k in range(3):
    M[k] = tmp[k]*np.sqrt(2*(E+g0*np.dot(r, gamma))/(np.dot(tmp, omega)))

print('Итого точка в координатах Кузнецова', M, gamma)
print(*calcIntegrals(M, g0, gamma))
