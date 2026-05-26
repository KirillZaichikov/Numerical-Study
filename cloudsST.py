'''
Построение облаков непрерывности, подсчет минимального угла, фазовый портрет, градиент угла от точки на аттракторе
'''

import numpy as np
import time
import math as m
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from params.params import *
from utils.linal import *
from utils.integrator import *

import pyvista as pv

np.set_printoptions(suppress=True)

def fix_error(der, points, skip, i):
    # print("ERROR!!!")
    # print("ERROR_SKIP_aftermap", oldold)
    # print("ERROR_SKIP_beforemod", der)
    if points[skip - 1 - i][2] > 0 and der[2] < 0:
        der[2] = der[2] + np.pi * 2
    elif points[skip - 1 - i][2] < 0 and der[2] > 0:
        der[2] = der[2] - np.pi * 2
    if np.linalg.norm(der) > 0.1:
        print('''LARGE_ERROR
                 LARGE_ERROR
                 LARGE_ERROR
                LARGE_ERROR
                LARGE_ERROR''')
    print("ERROR", i, der, points[skip - 1 - i])
    # print(np.linalg.norm(der))
    # print("old_der", old)
    # print("old_point", points[skip - i])

def func_for_get_to_attractor(params, start_point, skip_num, skip_var_num):
    start = time.time()
    for i in range(skip_num):
        makeStep(start_point, dimension, diffFunc, params, step)
    der = start_point.copy()
    der[0] += eps
    for i in range(skip_var_num):
        makeStep(start_point, dimension, diffFunc, params, step)
        makeStep(der, dimension, diffFunc, params, step)

        der = der - start_point
        new_norm = np.linalg.norm(der)
        der = (der / new_norm) * eps
        der = der + start_point

    last_vector_for_cu, last_vector_for_lyap = der, der
    end = time.time()
    print(f'time for get attractor: {end-start}')
    print("last_point", start_point)
    print("last_vector_for_cu", last_vector_for_cu)
    return start_point, last_vector_for_cu, last_vector_for_lyap

def find_min_ang(params, start_point, first_vector_for_cu, first_vector_for_lyap, count_of_trajectory, skip_var_num, points_main, points_skip, vectors_cu, vectors_ss):
    start = time.time()
    angles, min_ang = [None] * count_of_trajectory, 10
    der = first_vector_for_cu
    p1 = 0

    # ВЫЧИСЛЕНИЕ В ПРЯМОМ ВРЕМЕНИ
    for i in range(count_of_trajectory):
        makeStep(der, dimension, diffFunc, params, step)
        makeStep(start_point, dimension, diffFunc, params, step)

        points_main[i] = start_point.copy()

        der = der - start_point
        # print(np.linalg.norm(der))
        # КОСТЫЛЬ ДЛЯ РАСЧЕТОВ НА ТОРЕ
        if np.linalg.norm(der) > np.pi:
            if points_main[i][2] > 0 and der[2] < 0:
                der[2] = der[2] + np.pi * 2
            elif points_main[i][2] < 0 and der[2] > 0:
                der[2] = der[2] - np.pi * 2
            print("ERROR_DIRECT", i, der, points_main[i])
            # print(start_point)
            # print("old_der", old)
            # print("old_point", points_main[i][i-1])
    # if i > 15 and i < 20:
        #     print("der", der)
            # print("norm", np.linalg.norm(der))
        new_norm = np.linalg.norm(der)
        vectors_cu[i] = (der / new_norm).copy()
        der = (der / new_norm) * eps
        der = der + start_point

        # if i > 15 and i < 20:
        #     print("point", points_main[i])
        #     print("vector", vectors_cu[i])
        #     print("new_norm", new_norm)
        p1 += m.log(new_norm / eps)

    print('direct lyap', p1 / integrate_time)
    last_lyap_vector = der
    last_point = start_point
    last_vector = der

    # ЗАПОМИНАЕМ ТОЧКИ ДЛЯ ОБРАТНОГО
    for i in range(skip_var_num):
        points_skip[i] = start_point.copy()
        makeStep(start_point, dimension, diffFunc, params, step)

    der = np.identity(dimension) * eps
    for i in range(3):
        der[i] += start_point

    # ПРОГРЕВ В ОБРАТНОМ ВРЕМЕНИ
    for i in range(skip_var_num-1):
        for j in range(3):
            # old = der[j].copy()
            # print("before",der[j])
            makeStep(der[j], dimension, diffFuncRev, [*params], step)
            # makeStep(der[j], dimension, diffFuncRev, [*params, points_skip[skip_var_num - 1 - i][0]], step)
            # print(points_skip[skip_var_num - 1 - i])
            # oldold = der[j].copy()
            der[j] = der[j] - points_skip[skip_var_num - 1 - i]
            # print("after",der[j])
            if np.linalg.norm(der[j]) > np.pi:
                print("Skip", i, j)
                # print("old_copy", old)
                print(points_skip[skip_var_num - 1 - i])
                fix_error(der[j], points_skip, skip_var_num, i)

        ortVecs(der, dimension, lyap_num)

        for j in range(3):
            new_norm = np.linalg.norm(der[j])
            der[j] = (der[j] / new_norm) * eps
            der[j] = der[j] + points_skip[skip_var_num - 1 - i]

    # РАСЧЕТЫ В ОБРАТНОМ ВРЕМЕНИ
    p2 = [0,0,0]
    for i in range(count_of_trajectory):
        for j in range(3):
            # old = der[j].copy()
            makeStep(der[j], dimension, diffFuncRev, [*params], step)
            # makeStep(der[j], dimension, diffFuncRev, [*params, points_main[count_of_trajectory - 1 - i][0]], step)
            # print(f"vector{j}:",der[j])
            der[j] = der[j] - points_main[count_of_trajectory - 1 - i]
            if np.linalg.norm(der[j]) > np.pi:
                print(points_main[count_of_trajectory - 1 - i])
                print("CALC", i, j)
                fix_error(der[j], points_main, count_of_trajectory, i)
        # print(f"main_traj:",points_main[count_of_trajectory - 1 - i])

        ortVecs(der, dimension, lyap_num)
        for j in range(3):
            new_norm = np.linalg.norm(der[j])
            p2[j] += m.log(new_norm / eps)

        for j in range(3):
            der[j] = der[j] / np.linalg.norm(der[j])

        vectors_ss[count_of_trajectory - 1 - i] = der[2]
        angles[i] = abs(np.pi / 2 - np.arccos(
            np.dot(vectors_cu[count_of_trajectory - 1 - i], vectors_ss[count_of_trajectory - 1 - i])))
        if angles[i] < 0.001:
            print("ZERO ANGLE!!!!!", i)
        for j in range(3):
            der[j] = (der[j] / np.linalg.norm(der[j])) * eps
            der[j] = der[j] + points_main[count_of_trajectory - 1 - i]


    print('back time lyap 1:', p2[0] / integrate_time)
    print('back time lyap 2:', p2[1] / integrate_time)
    print('back time lyap 3:', p2[2] / integrate_time)

    end = time.time()
    print(f'time for num: {end-start}')
    min_ang = np.array(angles).min()
    return angles, min_ang, last_point, last_vector, last_lyap_vector

def draw_continuity_cloud(axCU, axSS, skip, points, vectors_ss, vectors_cu, points_len):
    angles, dist = [], []
    # vectors_ss = np.flip(vectors_ss, 0)
    print(points_len/skip)
    point1 = []
    point2 = []
    for i in range(int(points_len / skip)):
        for j in range(int(points_len / skip)):
            if j <= i:
                continue
            # angle = np.arccos(np.dot(vectors_ss[i*skip], vectors_ss[j*skip])/(np.linalg.norm(vectors_ss[i*skip])*np.linalg.norm(vectors_ss[j*skip])))
            angle = np.arccos(np.dot(vectors_ss[i*skip], vectors_ss[j*skip]))
            # angles.append(angle)
            if not np.isnan(angle):
                angles.append(angle)
            else:
                angles.append(0)
            dist.append(m.dist(points[i*skip], points[j*skip]))
            # if angles[-1] > 0.2 and angles[-1] < np.pi - 0.2 and dist[-1]<0.05:
            #     print(i, j)
            #     print(vectors_ss[i*skip])
            #     print(vectors_ss[j*skip])
            #     print(angles[-1])
            #     print(dist[-1])
            # if np.isnan(angles[-1]):
            #     print("ang", angles[-1])
            #     print("dist", dist[-1])
            #     print(vectors_ss[i*skip])
            #     print(vectors_ss[j*skip])
            #     point1.append(points[i*skip])
            #     point2.append(points[j*skip])
    
    point1_ = np.array(point1).T
    point2_ = np.array(point2).T
    print(point1_)
    print(len(dist), len(angles))
    axSS.scatter(dist, angles, c='#1F77B4', s=0.05)

    angles = []
    point1 = []
    point2 = []
    for i in range(int(points_len / skip)):
        for j in range(int(points_len / skip)):
            if j <= i:
                continue
            angle = np.arccos(np.dot(vectors_cu[i*skip], vectors_cu[j*skip])/(np.linalg.norm(vectors_cu[i*skip])*np.linalg.norm(vectors_cu[j*skip])))
            # angles.append(angle)
            if not np.isnan(angle):
                angles.append(angle)
            else:
                angles.append(0)
            # if angles[-1] > 0.5 and angles[-1] < np.pi - 0.5:
            #     print(i, j)
            #     print(vectors_cu[i*skip])
            #     print(vectors_cu[j*skip])
            #     print(m.dist(points[i*skip], points[j*skip]))
                # point1.append(points[i*skip])
                # point2.append(points[j*skip])
            # dist.append(m.dist(points[i*skip], points[j*skip]))
    print(len(dist), len(angles))
    axCU.scatter(dist, angles, c='#1F77B4', s=0.05)

    axSS.plot([eps, eps], [0, np.pi], c="red")
    axSS.set_xlim(0, max(dist))
    axSS.set_ylim(0, np.pi)
    axCU.plot([eps, eps], [0, np.pi], c="red")
    axCU.set_xlim(0, max(dist))
    axCU.set_ylim(0, np.pi)


num = 1
min_ang = 10
all_angles = []
last_point, last_cu_vector, last_vector_for_lyap = func_for_get_to_attractor(params, initial_point, skip_time, int(skip_var_time / step))
points_main = np.zeros(shape=(int(integrate_time / step / num), dimension), dtype=np.float32)
points_skip = np.zeros(shape=(int(skip_var_time / step), dimension), dtype=np.float32)
vectors_cu, vectors_ss = np.zeros(shape=(int(integrate_time / step / num), dimension), dtype=np.float32), np.zeros(shape=(int(integrate_time / step / num), dimension), dtype=np.float32)

for i in range(int(num)):
    angles, angle, last_point, last_cu_vector, last_vector_for_lyap = find_min_ang(params, last_point, last_cu_vector, last_vector_for_lyap,
                                                        int(integrate_time / step / num), int(skip_var_time/step), points_main, points_skip, vectors_cu, vectors_ss)
    if min_ang > angle:
        min_ang = angle
    all_angles.append(angles)

fig, axCu = plt.subplots()
fig1, axSs = plt.subplots()
fig.suptitle( 'U', fontsize=20)
fig1.suptitle('S', fontsize=20)
if with_clouds:
    draw_continuity_cloud(axCu, axSs, skip_for_phase, points_main, vectors_ss, vectors_cu, int(integrate_time / step / num))
print(f'min angle: {min_ang}')

####################### проверка битых точек
# mask = np.abs(vectors_ss[:, 0] - 1) > 0.0008
# P_bad = (points_main[mask]).T
# # print(P_bad)
# axPP.plot(P_bad[0], P_bad[1], P_bad[2], color = "red", linestyle="", marker="o", markersize=1) # отрисовка фазового портрета

####################### ОБЫЧНЫЙ ФАЗОВЫЙ
# figPP = plt.figure()
# axPP = figPP.add_subplot(projection='3d')
pT=points_main.T
# axPP.plot(pT[0], pT[1], pT[2], linestyle="", marker="o", markersize=0.6) # отрисовка фазового портрета
######################

# проверка точек с нулевым углом
npang = np.array(angles)[::-1]
print(npang.shape)
mask = npang < 0.001
zero_points = (points_main[mask]).T
if len(zero_points[0])>0:
    print("кол-во точек с нулевым углом", len(zero_points[0]))
    print("нулевые точки\n")
    print("ksi", zero_points[0])
    print("eta", zero_points[1])
    print("teta", zero_points[2])
    # axPP.plot(zero_points[0], zero_points[1], zero_points[2], color = "red", linestyle="", marker="o", markersize=4)

##################################### РАЗВЕРНУТЫЙ ФАЗОВЫЙ ПОРТРЕТ
fig = plt.figure(figsize=(10, 6))

# Создаем сетку 2×2
gs = fig.add_gridspec(2, 2)

# 3D-график (занимает всю левую колонку)
ax3d = fig.add_subplot(gs[:, 0], projection='3d')

# 2D-графики справа
ax2d_1 = fig.add_subplot(gs[0, 1])
ax2d_2 = fig.add_subplot(gs[1, 1])

# ax3d.set_title("3D поверхность sin(√(x² + y²))")
ax3d.set_xlabel("ksi")
ax3d.set_ylabel("eta")
ax3d.set_zlabel("teta")

# ax2d_1.set_title("ksi")
ax2d_1.set_xlabel("ksi")
ax2d_1.set_ylabel("eta")

# ax2d_2.set_title("2D scatter cos(x)")
ax2d_2.set_xlabel("teta")
ax2d_2.set_ylabel("eta")

every_iter = 10
ax3d.plot(pT[0][::every_iter], pT[1][::every_iter], pT[2][::every_iter], linestyle="", marker="o", markersize=0.6) # отрисовка фазового портрета
ax3d.plot(zero_points[0], zero_points[1], zero_points[2], color = "red", linestyle="", marker="o", markersize=3)

ax2d_1.plot(pT[0][::every_iter], pT[1][::every_iter], linestyle="", marker="o", markersize=0.3, color="black", rasterized=True)
ax2d_2.plot(pT[2][::every_iter], pT[1][::every_iter], linestyle="", marker="o", markersize=0.3, color="black", rasterized=True)

ax2d_1.plot(zero_points[0], zero_points[1], color = "red", linestyle="", marker="o", markersize=3)
ax2d_2.plot(zero_points[2], zero_points[1], color = "red", linestyle="", marker="o", markersize=3)

##################################################
#УГОЛ ГРАДИЕНТОМ

x_n = (pT[0] - pT[0].min()) / (pT[0].max() - pT[0].min())
y_n = (pT[1] - pT[1].min()) / (pT[1].max() - pT[1].min())
z_n = (pT[2] - pT[2].min()) / (pT[2].max() - pT[2].min())

points = np.column_stack((x_n, y_n, z_n))
cloud = pv.PolyData(points)
cloud["angle"] = npang

cmap = LinearSegmentedColormap.from_list(
    "yellow_red",
    ["#ffff00", "#ff0000"]
)

threshold = 0.001

pl = pv.Plotter()

# pl.show_axes()
pl.hide_axes()
pl.remove_bounds_axes()
# pl.show_bounds(
#     grid='front',
#     location='outer',
#     all_edges=True
# )

# --- точки с param > 0.001 (градиент) ---
pl.add_points(
    cloud.threshold(threshold, scalars="angle"),
    # scalars="angle",
    cmap=cmap,           # красный → жёлтый
    point_size=2,
    render_points_as_spheres=False
)

# --- точки с param <= 0.001 (синие) ---
if len(zero_points[0]) !=0:
    pl.add_points(
        cloud.threshold(threshold, invert=True, scalars="angle"),
        color="blue",
        point_size=5,
        render_points_as_spheres=False
    )

# pl.show_bounds(
#     bounds=(0, 1, 0, 1, 0, 1),
#     # grid='front',
#     # location='outer',
#     all_edges=True
# )
# pl.view_isometric()
# pl.show(auto_close=False, interactive_update=True)
# pl.show(auto_close=False)
# pl.add_scalar_bar(
#     title="Angle",
#     vertical=True,
#     position_x=0.85,
#     position_y=0.1,
#     width=0.08,
#     height=0.6,
#     label_font_size=12,
#     title_font_size=14
# )
pl.remove_scalar_bar()

def save_on_close():
    pl.screenshot("figure.png", window_size = [2000,2000])

pl.add_key_event("q", lambda:(save_on_close(),pl.close()))
pl.show(screenshot="figure")
# time.sleep(10)
print(pl.camera_position)
# pl.screenshot("figure.png", window_size = [2000,2000])
# mask_less = mask
# mask_more = npang > 0.001

# figAng = plt.figure()
# ax_ang = figAng.add_subplot(projection='3d')

# cmap = colors.LinearSegmentedColormap.from_list(
#     'yellow_red', ['yellow', 'red']
# )

# norm = colors.Normalize(
#     vmin=npang[mask_more].min(),
#     vmax=npang[mask_more].max()
# )

# sc = ax_ang.scatter(
#     pT[0][mask_more],
#     pT[1][mask_more],
#     pT[2][mask_more],
#     c=npang[mask_more],
#     cmap=cmap,
#     norm=norm,
#     s=1,
#     linewidths=0,
#     alpha=1
# )
# ax_ang.scatter(
#     pT[0][mask_less],
#     pT[1][mask_less],
#     pT[2][mask_less],
#     color='blue',
#     s=5,
#     linewidths=0,
#     alpha=1
# )

# cbar = plt.colorbar(sc, ax=ax_ang, pad=0.1)
# cbar.set_label('Параметр')
###################################################

###### СТРОЮ периодические точки
###### [0.05, 0.43, 0.2, 1]
# initial_guesses = [
#         # [np.array([-0.0316478,  -0.06491352, -2.92705166]), 1],
#         [np.array([0.03395866, 0.01807846, 0.87946063]), 5],
#         [np.array([0.00893707, -0.05983158, -2.62629825]), 7],
#         [np.array([0.0243333,   0.01052277, -1.84813452]), 6],
#         [np.array([-0.00345698,  0.06176062,  2.1358999 ]), 3],
#         [np.array([0.00174419, 0.05939727, 2.15595705]), 5],
#         [np.array([0.02218712,  0.00897619, -1.8092251]), 4],
#         [np.array([-0.07020054, -0.03911589, -0.30096515]), 5],
#         [np.array([-0.01830785, -0.06929036,  0.18548144]), 6]
# ]
# initial_guesses = [[np.array([0.02693017, -0.02482757, -2.59385737]), 2],
#                    [np.array([0.04050016, -0.01855747, -2.32504418]), 4],
#                    [np.array([0.06967844, -0.02299013, -2.30025628]), 7]]

# periodic_points = []
# for i in range(len(initial_guesses)):
#     x = initial_guesses[i][0].copy()
#     periodic_points.append([])
#     for j in range(initial_guesses[i][1]):
#         makeStep(x, dimension, diffFunc, params, step)
#         periodic_points[i].append(x.copy())

# for a in periodic_points:
#     coord = np.array(a).T
#     print(coord)
#     ax3d.plot(coord[0], coord[1], coord[2], linestyle="", marker="o", color='black', markersize=3)

# неподвижные точки
###### [0.05, 0.43, 0.2, 1]
# ax3d.plot(-0.02705276, -0.06695892, 0.17971049, linestyle="", marker="o", color='magenta', markersize=3)
# ax3d.plot(-0.0316478, -0.06491352, -2.92705166, linestyle="", marker="o", color='magenta', markersize=3)

plt.show()