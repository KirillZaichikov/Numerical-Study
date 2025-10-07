import math as m
import numpy as np

# M = np.ndarray(3)
# gamma = np.ndarray(3)
M = np.array([176.39671645238903, -168.85257112196294, -94.877918535603939])
gamma = np.array([0.69524430378189617, 0.40429993079272591, -0.59428864887075838])
gamma = -gamma

L = M[2]
G = m.sqrt(M[0]**2+M[1]**2+M[2]**2)
H = np.dot(M,gamma)
l = np.atan(M[1]/M[0])
g = np.atan((M[1]*gamma[0]-M[0]*gamma[1])/(H*L/G-G*gamma[2]))

print(L,G,H,l,g, sep='\n')

l = 3.6547925893
LG = 0.6508187558
HG = -0.3626716311

gamma[0] = (HG*m.sqrt(1-LG**2)+LG*m.sqrt(1-HG**2)*m.cos(0))*m.sin(l)+m.sqrt(1-HG**2)*m.sin(0)*m.cos(l)
gamma[1] = (HG*m.sqrt(1-LG**2)+LG*m.sqrt(1-HG**2)*m.cos(0))*m.cos(l)-m.sqrt(1-HG**2)*m.sin(0)*m.cos(l)
gamma[2] = LG*HG-m.sqrt(1-LG**2)*m.sqrt(1-HG**2)*m.cos(0)