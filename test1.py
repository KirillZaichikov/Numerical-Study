import matplotlib.pyplot as plt
import numpy as np


# mu = 0.104
# mu = 0.07157157157157157 # первая гетероклиника
# mu = 0.02702702702702703 # вторая гетероклиника
# mu = 0.012262262262262236 # третья гетероклиника
mu = 0.006506506506506526 # 4 гетероклиника
A = 0.95
# nu = 0.964
# nu = 0.8568568568568569 # первая гетероклиника
# nu = 0.9459459459459459 # вторая гетероклиника
# nu = 0.9754754754754755 # третья гетероклиника
nu = 0.986986986986987 # 4 гетероклиника


mu = 0.17393320964749534
nu = 0.8758103537981272
LE_1 = 0.2215690077628535
time = 4
times = 1
skip_time = 0
x_ = 0
lx, ly = [], []

def map(x):
    return -mu+A*abs(x)**nu

def map_rev(x):
    return abs(((x+mu)/A)**(1/nu))


for i in range(time):
    lx.append(x_)
    for i in range(times):
        x_ = map(x_)
    ly.append(x_)


fig, ax = plt.subplots()
mas_fx = np.linspace(min(lx[skip_time:]), max(lx[skip_time:]), 3000)

mas_fy = -mu+A*abs(mas_fx)**nu
# ax.plot(mas_fx, mas_fy, c='magenta')

# mas_fy = -mu+A*abs(-mu+A*abs(mas_fx)**nu)**nu
# ax.plot(mas_fx, mas_fy, c='green')

# l1 = np.array(lx)
# l2 = np.array(ly)

ax.scatter(lx[skip_time:],ly[skip_time:], s=0.1, c='black') # s = 0.1
ax.scatter(lx[-1],ly[-1], s=5, c='red') # s = 0.1
ax.plot(mas_fx,mas_fx, c='red')
plt.show()