import matplotlib.pyplot as plt
import numpy as np
import time as t
import math as m


A = 0.95
time = 3
skip_time = 0
lx, ly = [], []

def map(x, nu, mu):
    return -mu+A*abs(x)**nu

angle = np.linspace(-0.2, -5, 500)
critical_points_nu = np.zeros((len(angle), 10))
critical_points_mu = np.zeros((len(angle), 10))
for line_num, a in enumerate(angle):
    count = 0 
    times = 1
    print(1-m.sqrt(0.08/(1+a**2)))
    x = np.linspace(1-m.sqrt(0.05/(1+a**2)), 1, 1000)
    y = (x-1)*(a)
    neg = True
    for index, element in enumerate(x):
        x_ = 0.0
        for i in range(time):
            for j in range(times):
                x_ = map(x_, element, y[index])

        x_tmp = x_
        for j in range(times):
            x_tmp = map(x_tmp, element, y[index])

        if neg:
            if x_ <  x_tmp: 
                times = times * 2
                neg=False
                if times > 1024:
                    break
                critical_points_nu[line_num][count] = element
                critical_points_mu[line_num][count] = y[index]
                count += 1
        else:
            if x_ >  x_tmp: 
                times = times * 2
                neg=True
                if times > 1024:
                    break
                critical_points_nu[line_num][count] = element
                critical_points_mu[line_num][count] = y[index]
                count += 1


np.savetxt('mas.txt', critical_points_nu.T[:5])
fig, ax = plt.subplots()
# mas = np.array(critical_points).T
ax.plot(critical_points_nu.T[:3].T, critical_points_mu.T[:3].T, c='black')
ax.plot([1, 1], [0, 1], c='magenta')
plt.show()
