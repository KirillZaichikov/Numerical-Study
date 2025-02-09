import numpy as np
import matplotlib.pyplot as plt


def Shimizu_3D_flow(state, res, params, H) -> None:
    res[0] = (state[1]) / H
    res[1] = (- params[0] * state[1] - state[0] * state[2] + state[0]) / H
    res[2] = (- params[1] * state[2] + pow( state[0] , 2.0 )) / H

def ShimizuX3_3D_flow(state, res, params, H) -> None:
    # Param - alpha lamda B 
    res[0] = state[1] / H
    res[1] = (params[2] * pow ( state[0] , 3.0 ) - params[0] * state[1] - state[0] * state[2] + state[0]) / H
    res[2] = (- params[1] * state[2] + pow ( state[0] , 2.0 )) / H

def dverkStep(val, dimension, diffFunc, params, step, H=1) -> None:
    k1, k2, k3, k4, k5, k6, k7, k8 = np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), \
        np.zeros(dimension), np.zeros(dimension), np.zeros(dimension), np.zeros(dimension)
    arg = np.zeros(dimension)
    diffFunc(val, k1, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step / 6. * k1[j]
    
    diffFunc(arg, k2, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j])
    
    diffFunc(arg, k3, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j])
    
    diffFunc(arg, k4, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j])
    
    diffFunc(arg, k5, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j])
    
    diffFunc(arg, k6, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j])
    
    diffFunc(arg, k7, params, H)
    for j in range(dimension):
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j])
    
    diffFunc(arg, k8, params, H)
    for j in range(dimension):
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j])


alpha = 1.0470925324675324
lamda = 0.5988294354838709
LE_1 = 0.0026132984145119
crossesction = 1
params = np.array(1.0470925324675324, 0.5988294354838709, 0.2]) # alpha lamda B
# params = np.array([1.13, 0.48]) # alpha lamda 
initial_point = np.array([0.000000001, 0, 0])
time = 5000
step = 0.01
dimension = 3
coord_list = []
coord_list_for_phase = []

print('1) saddle index is ', -((-params[0]+(params[0]**2+4)**(1/2))/2) / (-params[1]))
nu = -((-params[1]) / ((-params[0]+(params[0]**2+4)**(1/2))/2))
print('2) saddle index is ', nu)

f = False
# итерации потока
for i in range(int(time/step)):
    z_last = initial_point[2] - crossesction
    dverkStep(initial_point, dimension, ShimizuX3_3D_flow, params, step)
    if (initial_point[2] - crossesction > 0) and (z_last < 0) and not f:
        H = - params[1] * initial_point[2] + initial_point[0] ** 2
        dverkStep(initial_point, dimension, ShimizuX3_3D_flow, params, -(initial_point[2] - crossesction), H)
        coord_list.append(list(initial_point[:3]))
        f = True
    elif (initial_point[2] - crossesction > 0) and (z_last < 0) and f:
        f=False
    coord_list_for_phase.append(list(initial_point))

# отрисовка фазового
transpose_list1 = (np.array(coord_list_for_phase)).T
transpose_list = (np.array(coord_list)).T
fig1 = plt.figure()
ax1 = fig1.add_subplot(projection='3d')
ax1.plot(transpose_list1[0], transpose_list1[1], transpose_list1[2])
ax1.scatter(transpose_list[0], transpose_list[1], transpose_list[2], s=5, c="red")

# отрисовка сечения
fig, ax = plt.subplots()
ax.scatter(transpose_list[0], transpose_list[1], s=5, c="red")

# подготовка массивов для одномерного отображения
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

print('mu (max) ', max_point_y)
print('x0 (shift) ', max_point_x)
print('lenght of the poicare list ', len(lx))
print('angle is ', (max_point_y-min_point_y)/(max_point_x-min_point_x))

with open('firstretMaps\\'+'Ha_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:6]+'.txt', 'w') as dat:
    dat.write(str(nu)) # nu
    dat.write('\n')
    dat.write(str(max_point_x)) # shift
    dat.write('\n')
    dat.write(str(max_point_y)) # mu
    dat.write('\n')
    for x in lx:
        dat.write(str(x))
        dat.write(' ')
    dat.write('\n')
    for y in ly:
        dat.write(str(y))
        dat.write(' ')

plt.show()
