import math as m
import numpy as np

# M = np.ndarray(3)
# gamma = np.ndarray(3)
M = np.array([-54.065601,  -41.130199,  76.661598]) 
gamma = np.array([  0.267023,  0.203137,  0.942037])
gamma = -gamma # Нужно гамму обратно повернуть
d = 0.423

#Нужно обратно повернуть
QK = np.array([[np.cos(-d), np.sin(-d), 0],
                 [-np.sin(-d),np.cos(-d),0],
                 [0,0,1]])
M = QK @ M
gamma = QK @ gamma
# M = M / np.linalg.norm(M)

L = M[2]
G = m.sqrt(M[0]**2+M[1]**2+M[2]**2)
H = np.dot(M,gamma)
l = m.atan(M[0]/M[1])
g = np.atan((M[1]*gamma[0]-M[0]*gamma[1])/(H*L/G-G*gamma[2]))

print(L,G,H,l,g, sep='\n')
print("L/G", L/G, "H/G", H/G, "l", l+np.pi, sep='\n')
