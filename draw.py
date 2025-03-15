import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
params = np.array([1.0501456818181818, 0.6030381935483871, 0.2], dtype=np.longdouble) # alpha lamda B
m = np.load('2dimmaps/'+'a_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:9]+'.npy')

m = m.T
r = []
for a in m:
    if a[1]>0.4848:
        r.append(a)
r1 = np.array(r)
print(r1)
# r1 = r1.T
coord_list = r1

lx, ly = [], []
max_point_x, max_point_y = -100, -100
min_point_x, min_point_y = -100, -100
for i in range(len(coord_list)-1):
    if coord_list[i][0] > 0.0 and (abs(coord_list[i+1][0]) > coord_list[i][0]): # второе условие нужно если есть лакуна
        lx.append(coord_list[i][0])
        ly.append(abs(coord_list[i+1][0])) # abs
        if abs(coord_list[i+1][0]) > max_point_y:
            max_point_x = coord_list[i][0]
            max_point_y = abs(coord_list[i+1][0])
        elif coord_list[i][0] > min_point_x:
            min_point_x = coord_list[i][0]
            min_point_y = abs(coord_list[i+1][0])

print(lx[0])
ax.scatter(lx, ly, s=1, c='red')
plt.show()
