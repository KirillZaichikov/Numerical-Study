import math as m
import numpy as np

# M = np.ndarray(3)
# gamma = np.ndarray(3)
M = np.array([-62.6672, -30.7099, 57.8557]) 
gamma = np.array([0.269143, 0.131893, 0.954026]) #- Это неподвижная точка 0.485 752
gamma = -gamma # Нужно гамму обратно повернуть
M = np.array([-54.06849765, -41.5933255,   76.16298962]) 
gamma = np.array([0.26638367, 0.20492122, 0.94183174])
gamma = -gamma
d = 0.485

#Нужно обратно повернуть
QK = np.array([[np.cos(-d), np.sin(-d), 0],
                 [-np.sin(-d),np.cos(-d),0],
                 [0,0,1]])
M = QK @ M
gamma = QK @ gamma
print(M,gamma)
# M = M / np.linalg.norm(M)

L = M[2]
G = m.sqrt(M[0]**2+M[1]**2+M[2]**2)
H = np.dot(M,gamma)
l = m.atan(M[0]/M[1])
g = np.atan((M[1]*gamma[0]-M[0]*gamma[1])/(H*L/G-G*gamma[2]))

print(L,G,H,l,g, sep='\n')
print("L/G", L/G, "H/G", H/G, "l", l+np.pi, sep='\n')
