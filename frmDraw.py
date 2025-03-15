import numpy as np
import matplotlib.pyplot as plt


lx_tmp = []
ly_tmp = []
params = np.array([1.371672077922078, 0.43540322580645163, 0.2], dtype=np.longdouble) # alpha lamda B   HIGH PRIORITY 0
params = np.array([1.05167, 0.6013812096774193, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 2
params = np.array([1.364,0.459, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ПЕРВОГО КРИТЕРИЯ Лакуна 2
# params = np.array([1.045,0.588, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ВТОРОГО КРИТЕРИЯ Лакуна 4
# params = np.array([1.04476, 0.60583, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 3
with open('firstretMaps\\'+'Ha_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:9]+'.txt', 'r') as dat:
# with open('firstretMaps\\HighPriority\\'+'Ha_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:9]+'.txt', 'r') as dat:
    nu = np.float64(dat.readline()) # separatrix value
    shift = np.float64(dat.readline()) # shift
    mu = np.float64(dat.readline()) # mu
    shift_min = np.float64(dat.readline()) # for cline is down
    mu_min = np.float64(dat.readline()) # for cline is down
    tmp = dat.readline()
    tmp = (tmp.split(' '))[:-1]
    for x in tmp:
        lx_tmp.append(np.float64(x))
    tmp = dat.readline()
    tmp = (tmp.split(' '))[:-1]
    for y in tmp:
        ly_tmp.append(np.float64(y))

print('shift', shift)
print('nu', nu)
print('mu', mu)
print('shift_min', shift_min)
print('mu_min', mu_min)

fig, ax = plt.subplots(figsize=(5, 5))
lx = np.array(lx_tmp)
ly = np.array(ly_tmp)
print(len(lx))
ax.scatter(lx, ly, s=2, c='black')

x = np.linspace(min(lx), max(lx), 3000)
y = abs(0.51*pow((abs(x-shift)), nu)-mu)
# y = abs(1.35*pow((abs(x-shift_min)), nu)+mu_min) # for cline is down
ax.plot(x, y, c='red')
ax.plot(x, x, c='blue', linestyle='--')

plt.show()
