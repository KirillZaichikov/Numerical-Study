import matplotlib.pyplot as plt
import numpy as np


mu = 1#0.104
A = 0.2
C = 1.5
nu = 0.5#0.964
time = 10000
times = 1
skip_time = 0
x_ = 0.0
lx, ly = [], []

# def map(x):
#     return -mu+A*abs(x)**nu

def map(x):
    return abs(-mu+A*abs(x)**nu+C*(abs(x)**(2*nu)))

for i in range(time):
    lx.append(x_)
    for i in range(times):
        x_ = map(x_)
    ly.append(x_)


fig, ax = plt.subplots()
# mas_fx = np.linspace(min(lx[skip_time:]), max(lx[skip_time:]), 3000)

mas_fx = np.linspace(0., 0.8, 3000)

mas_fy = -mu+A*abs(mas_fx)**nu
# ax.plot(mas_fx, mas_fy, c='magenta')

mas_fy = -mu+A*abs(-mu+A*abs(mas_fx)**nu)**nu
# ax.plot(mas_fx, mas_fy, c='green')


mas_sec_ord = abs(-mu+A*abs(mas_fx)**nu+C*abs(mas_fx)**(2*nu))
ax.plot(mas_fx, mas_sec_ord, c='green')

mas_sec_ord = abs(-mu+A*(abs(-mu+A*abs(mas_fx)**nu+C*abs(mas_fx)**(2*nu)))**nu+C*abs((-mu+A*abs(mas_fx)**nu+C*abs(mas_fx)**(2*nu)))**(2*nu))
ax.plot(mas_fx, mas_sec_ord, c='magenta')
# l1 = np.array(lx)
# l2 = np.array(ly)

ax.scatter(lx[skip_time:],ly[skip_time:], s=4, c='black')
ax.plot(mas_fx,mas_fx, c='red')
plt.show()

