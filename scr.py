import matplotlib.pyplot as plt
import numpy as np

mu = 0.2#0.104
A = 0.5305
nu = 0.7#0.964
time = 20000
times = 1
skip_time = 0
x_ = 0.0
lx, ly = [], []

def map(x):
    return -mu+A*abs(x)**nu


for i in range(time):
    lx.append(x_)
    for i in range(times):
        x_ = map(x_)
    ly.append(x_)


fig, ax = plt.subplots()
mas_fx = np.linspace(min(lx[skip_time:]), max(lx[skip_time:]), 3000)

mas_fy = -mu+A*abs(mas_fx)**nu
ax.plot(mas_fx, mas_fy, c='magenta')

mas_fy = -mu+A*abs(-mu+A*abs(mas_fx)**nu)**nu
# ax.plot(mas_fx, mas_fy, c='green')

# l1 = np.array(lx)
# l2 = np.array(ly)

ax.scatter(lx[skip_time:],ly[skip_time:], s=0.1, c='black')
ax.plot(mas_fx,mas_fx, c='red')
plt.show()

