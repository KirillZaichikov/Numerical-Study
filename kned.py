import matplotlib.pyplot as plt
import numpy as np
import time
import distinctipy

mu = np.linspace(0, 0.5, 500)
A = 0.63
nu = np.linspace(0.5, 1, 500)
# times = 1
x_ = 0.0
lx, ly = [], []
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(nu[0], nu[-1])
ax.set_ylim(mu[0], mu[-1])

def map(x, mu, nu):
    return -mu+A*abs(x)**nu

lines = [''] * len(mu) * len(nu)
for i, par1 in enumerate(nu):
    for j, par2 in enumerate(mu):
        x_=0
        x_ = map(x_, par2, par1)
        while len(lines[i*len(nu)+j]) < 11:
            x_ = map(x_, par2, par1)
            if x_ > 0:
                lines[i*len(nu)+j] += '1'
            else:
                lines[i*len(nu)+j] += '0'

uniq = set(lines)
print(len(uniq))
uniq_l = list(uniq)
del(uniq)
start = time.time()
colors = distinctipy.get_colors(len(uniq_l))
print(time.time()-start)
colormap = np.array(colors)
x, y, categories= [], [], []

for i in range(len(nu)):
    for j in range(len(mu)):
        x.append(nu[i])
        y.append(mu[j])
        for a in uniq_l:
            if a == lines[i*len(nu)+j]:
                categories.append(uniq_l.index(a))
                break
ax.scatter(x, y, c=colormap[categories], s=1.2)
plt.show()