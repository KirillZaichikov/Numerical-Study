import matplotlib.pyplot as plt
import numpy as np
import time
import distinctipy
import colors.colors as clrs

DISCR = 1000
mu = np.linspace(0, 0.25, DISCR)
A = 0.5305
nu = np.linspace(0.5, 1, DISCR)
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
        while len(lines[i*len(nu)+j]) < 12:
            x_ = map(x_, par2, par1)
            if x_ > 0:
                lines[i*len(nu)+j] += '1'
            else:
                lines[i*len(nu)+j] += '0'

uniq = set(lines)
print("uniq seq: ", len(uniq))
uniq_l = sorted(list(uniq))
print(uniq_l[0])
del(uniq)
# start = time.time()
# colors = distinctipy.get_colors(len(uniq_l))
# print(time.time()-start)
colormap = np.array(clrs.colors800[:len(uniq_l)])
x, y, categories= [], [], []

for i in range(len(nu)):
    for j in range(len(mu)):
        x.append(nu[i])
        y.append(mu[j])
        for k, a in enumerate(uniq_l):
            if a == lines[i*len(nu)+j]:
                categories.append(uniq_l.index(a))
                break
            if k + 1 == len(uniq_l):
                categories.append(None)
ax.scatter(x, y, c=colormap[categories], s=0.1)
# plt.savefig('grafik.pdf', format='pdf')
plt.show()