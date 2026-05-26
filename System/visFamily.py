import numpy as np
from params import *
import matplotlib.pyplot as plt

def get_omega(start_point, params):
    ro, mpar, i1, i2, a0, g, eps = params
    dzeta, phi, teta = start_point

    omega2 = 2*(eps - mpar * g * np.sin(teta) * (ro + a0*np.sin(phi)) ) / \
        (i1 * np.cos(dzeta-phi) ** 2 + i2 * np.sin(dzeta-phi) ** 2 + \
            mpar * np.cos(teta)**2*(ro*np.cos(dzeta)-a0 * np.sin(dzeta-phi))**2)
    omega = m.sqrt(omega2)
    if omega2 < 0: print("ERROR: omega less than 0")
    return omega


phi = np.linspace(0, 2*np.pi, DISCR)
points = []
for _phi in phi:
    initial = np.array([0.01, _phi, 0.01])
    dzeta_, phi_, teta_ = initial
    omega = get_omega([0.01, _phi, 0.01], params)
    points.append(np.array([omega * np.cos(dzeta_-phi_), omega * np.sin(dzeta_-phi_),phi_]))

points = np.array(points).T
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.plot(points[2], points[0], points[1], linestyle="", marker="o", markersize=0.1, color="black")


plt.show()
