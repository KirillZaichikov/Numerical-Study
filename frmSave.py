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


crossesction = 1.25
# params = np.array([1.371672077922078, 0.43540322580645163, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 0
params = np.array([1.05167, 0.6013812096774193, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 1
params = np.array([1.04476, 0.60583, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 3
params = np.array([1.0498470941558442,0.6015436612903226, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 
params = np.array([1.04200000,0.59800000, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY
params = np.array([1.364,0.459, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ПЕРВОГО КРИТЕРИЯ Лакуна 2
params = np.array([1.045,0.588, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ВТОРОГО КРИТЕРИЯ Лакуна 4
params = np.array([1.049, 0.5961, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ВТОРОГО КРИТЕРИЯ Лакуна 8
params = np.array([1.06, 0.5925, 0.2], dtype=np.longdouble) # точка 2 Лакуна 8 A = 0.94
params = np.array([1.058, 0.594, 0.2], dtype=np.longdouble) # точка 2 Лакуна 16 A = 0.95
params = np.array([1.053, 0.599, 0.2], dtype=np.longdouble) # точка 2 Лакуна 32 A = 0.9
params = np.array([1.05, 0.6015, 0.2], dtype=np.longdouble) # точка 2 Лакуна 32 A = 0.978768
params = np.array([1.04475, 0.60575, 0.2], dtype=np.longdouble) # точка 3 Лакуна 256 A = 1.0176
params = np.array([1.04475, 0.6057, 0.2], dtype=np.longdouble) # точка 3 Лакуна 256 A = 1.002
params = np.array([1.351, 0.48, 0.2], dtype=np.longdouble) # alpha lamda ФИНАЛЬНАЯ ОЦЕНКА ВОЗЛЕ ПЕРВОГО КРИТЕРИЯ Лакуна 4
params = np.array([0.74221, 0.63016, 0.2], dtype=np.longdouble) # Счетное
params = np.array([0.697, 0.61, 0.2], dtype=np.longdouble) 
params = np.array([1.044, 0.605, 0.2], dtype=np.longdouble) # точка 3 Лакуна 64 A = 1.0176
iter_param = 4
initial_point = np.array([0.0000001, 0, 0], dtype=np.longdouble)
skip_time = 0
time = 400000
skip_phase = 5
step = 0.01
dimension = 3
coord_list_for_frm = []
coord_list_for_poincare = []
coord_list_for_phase = []

print('1) saddle index is ', -((-params[0]+(params[0]**2+4)**(1/2))/2) / (-params[1]))
nu = -((-params[1]) / ((-params[0]+(params[0]**2+4)**(1/2))/2))
print('2) saddle index is ', nu)

for i in range(int(skip_time/step)):
    dverkStep(initial_point, dimension, ShimizuX3_3D_flow, params, step)

# итерации потока
# iter_count = iter_param - 1  # старуем с iter_param - 1 чтобы взять именно боковые точки
iter_count = 0  # старуем с 0 чтобы иметь норм треугольник
for i in range(int(time/step)):
    z_last = initial_point[2] - crossesction
    dverkStep(initial_point, dimension, ShimizuX3_3D_flow, params, step)
    if (initial_point[2] - crossesction > 0) and (z_last < 0):
        H = - params[1] * initial_point[2] + initial_point[0] ** 2
        dverkStep(initial_point, dimension, ShimizuX3_3D_flow, params, -(initial_point[2] - crossesction), H)
        iter_count += 1
        if iter_count == iter_param:
            iter_count = 0
            coord_list_for_frm.append(list(initial_point[:3]))
        coord_list_for_poincare.append(list(initial_point[:2]))
    # if i % skip_phase == 0:
    #     coord_list_for_phase.append(list(initial_point))

# отрисовка фазового
transpose_list1 = (np.array(coord_list_for_phase, dtype=np.longdouble)).T
transpose_points_for_frm = (np.array(coord_list_for_frm, dtype=np.longdouble)).T
transpose_points_for_poincare = (np.array(coord_list_for_poincare, dtype=np.longdouble)).T
fig1 = plt.figure()
ax1 = fig1.add_subplot(projection='3d')
# ax1.plot(transpose_list1[0], transpose_list1[1], transpose_list1[2]) # отрисовка фазового портрета
ax1.scatter(transpose_points_for_frm[0], transpose_points_for_frm[1], transpose_points_for_frm[2], s=5, c="red") # отрисовка точек сечения которые идут для отображения первого возвращения
np.save('2dimmaps/'+'a_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:9], transpose_points_for_poincare)

# отрисовка сечения
fig, ax = plt.subplots()
ax.scatter(transpose_points_for_poincare[0], transpose_points_for_poincare[1], s=5, c="black") # отрисовка точек сечения
ax.scatter(transpose_points_for_frm[0], transpose_points_for_frm[1], s=5, c="red") # отрисовка точек сечения которые идут для отображения первого возвращения

# подготовка массивов для одномерного отображения
lx, ly = [], []
max_point_x, max_point_y = -100, -100
min_point_x, min_point_y = 1000, 1000
for i in range(len(coord_list_for_frm)-1): # пробую y вместо x
    # if (abs(coord_list_for_frm[i+1][0]) > coord_list_for_frm[i][0]): # условие нужно если есть лакуна
        lx.append(abs(coord_list_for_frm[i][0])) # здесь пробую модуль
        ly.append(abs(coord_list_for_frm[i+1][0])) # abs
        if abs(coord_list_for_frm[i+1][0]) > max_point_y:
            max_point_x = abs(coord_list_for_frm[i][0])
            max_point_y = abs(coord_list_for_frm[i+1][0])
        elif abs(coord_list_for_frm[i+1][0]) < min_point_y:
            min_point_x = abs(coord_list_for_frm[i][0])
            min_point_y = abs(coord_list_for_frm[i+1][0])

print('mu (max) ', max_point_y)
print('x0 (shift) ', max_point_x)
print('lenght of the poicare list ', transpose_points_for_poincare[0])
print('lenght of the акь list ', len(lx))
print('angle is ', (max_point_y-min_point_y)/(max_point_x-min_point_x))

fig2, ax2 = plt.subplots()
ax2.scatter(lx, ly, s=5, c="red")

with open('firstretMaps\\'+'Ha_'+f'{params[0]}'[:6]+'l_'+f'{params[1]}'[:9]+'.txt', 'w') as dat:
    dat.write(str(nu)) # nu
    dat.write('\n')
    dat.write(str(max_point_x)) # shift
    dat.write('\n')
    dat.write(str(max_point_y)) # mu
    dat.write('\n')
    dat.write(str(min_point_x)) # for cline down
    dat.write('\n')
    dat.write(str(min_point_y)) # for cline down
    dat.write('\n')
    for x in lx:
        dat.write(str(x))
        dat.write(' ')
    dat.write('\n')
    for y in ly:
        dat.write(str(y))
        dat.write(' ')

plt.show()
