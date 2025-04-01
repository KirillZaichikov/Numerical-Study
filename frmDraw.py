import numpy as np
import matplotlib.pyplot as plt


lx_tmp = []
ly_tmp = []
params = np.array([1.371672077922078, 0.43540322580645163, 0.2], dtype=np.longdouble) # alpha lamda B   HIGH PRIORITY 0
params = np.array([1.05167, 0.6013812096774193, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 2
params = np.array([1.364,0.459, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ПЕРВОГО КРИТЕРИЯ Лакуна 2
params = np.array([1.351,0.48, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ПЕРВОГО КРИТЕРИЯ Лакуна 4
params = np.array([1.395,0.415, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ПЕРВОГО КРИТЕРИЯ Лакуна 1
params = np.array([1.06, 0.5925, 0.2], dtype=np.longdouble) # точка 2 Лакуна 4 A = 0.94
params = np.array([1.058, 0.594, 0.2], dtype=np.longdouble) # точка 2 Лакуна 8 A = 0.95
params = np.array([1.056, 0.596, 0.2], dtype=np.longdouble) # точка 2 Лакуна 16 A = 0.9
params = np.array([1.053, 0.599, 0.2], dtype=np.longdouble) # точка 2 Лакуна 32 A = 0.9
params = np.array([1.0425, 0.604, 0.2], dtype=np.longdouble) # точка 3 Лакуна 32 A = 1.044
params = np.array([1.05, 0.6015, 0.2], dtype=np.longdouble) # точка 2 Лакуна 32 A = 0.978768
params = np.array([1.044, 0.605, 0.2], dtype=np.longdouble) # точка 3 Лакуна 256 A = 1.0176
# params = np.array([1.059, 0.593, 0.2], dtype=np.longdouble) 
# params = np.array([1.045,0.588, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ВТОРОГО КРИТЕРИЯ Лакуна 4
# params = np.array([1.049, 0.5961, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ВТОРОГО КРИТЕРИЯ Лакуна 8
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

x = np.linspace(min(lx), max(lx), 300000)
y = abs(1.0176*pow((abs(x-shift)), nu)-mu)
# y = abs(0.978768*pow((abs(x-shift_min)), nu)+mu_min) # for cline is down
ax.plot(x, y, c='red')
ax.plot(x, x, c='blue', linestyle='--')

plt.show()
