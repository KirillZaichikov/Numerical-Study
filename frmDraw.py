import numpy as np
import matplotlib.pyplot as plt


lx_tmp = []
ly_tmp = []
# params = np.array([1.372380481961717, 0.38544840620951815, 0.2]) # alpha lamda B
# params = np.array([1.0327122652063017, 0.5837562912782179, 0.2]) # alpha lamda B
# params = np.array([1.0429612661391086, 0.5976741914247489, 0.2]) # alpha lamda B
# params = np.array([1.0474084489589393,0.6027620038650215, 0.2]) # alpha lamda B
# params = np.array([1.0211057513914656,0.6093805069124424, 0.2]) # alpha lamda B
# params = np.array([1.021258116883117, 0.6092816935483871, 0.2]) # alpha lamda B
with open('firstretMaps\\'+'a_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:6]+'.txt', 'r') as dat:
    nu = np.float64(dat.readline()) # separatrix value
    shift = np.float64(dat.readline()) # shift
    mu = np.float64(dat.readline()) # mu
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

fig, ax = plt.subplots()
lx = np.array(lx_tmp)
ly = np.array(ly_tmp)
ax.scatter(lx, ly, s=4, c='black')

x = np.linspace(min(lx), max(lx), 3000)
y = abs(1.08*pow((abs(x-shift)), nu)-mu)
ax.plot(x, y)

plt.show()
