import matplotlib.pyplot as plt
from utils.integrator import dverkStep
import numpy as np
import distinctipy

z = np.linspace(0, 0.2, 100)
wl = []
for i, zp in enumerate(z):
    eps = -2 * zp
    # print(i)
    if len(z)/2:
        step = ((eps) * 2)/25
        colvo = 25
    else:
        step = ((eps) * 2)/50
        colvo = 50
    print(eps, step)
    wl.append([])
    for j in range(colvo):
        wl[i].append(-eps + step*j)

fig, ax = plt.subplots()

for xl in wl:
    for x in xl:
        ax.plot(xl, [z[wl.index(xl)]]*len(xl), c='black')

plt.show()