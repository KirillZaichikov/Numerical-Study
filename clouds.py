import numpy as np
import time
import math as m
import matplotlib.pyplot as plt
from params.params import *
from utils.linal import *
from utils.integrator import *




def func_for_get_to_attractor(params, start_point, skip_num, skip_var_num):
    start = time.time()
    der, der1 = [eps, 0, 0], [eps, 0, 0]
    der_norm, der_norm1 = der / np.linalg.norm(der), der1 / np.linalg.norm(der1)
    print(type(der_norm))
    for i in range(skip_num):
        makeStep(start_point, dimension, ShimizuX3_3D_flow, params, step)
    for i in range(skip_var_num):
        makeStepVar(der_norm, dimension, ShimizuX3_3D_flow_var, params, step, start_point)
        makeStepVar(der_norm1, dimension, ShimizuX3_3D_flow_var_trans_rev, params, step, start_point)
        makeStep(start_point, dimension, ShimizuX3_3D_flow, params, step)
        new_norm = np.linalg.norm(der_norm)
        new_norm1 = np.linalg.norm(der_norm1)
        der_norm = der_norm / new_norm
        der_norm1 = der_norm1 / new_norm1

    last_vector_for_cu, last_vector_for_lyap = der_norm1, der_norm
    end = time.time()
    print(f'time for get attractor: {end-start}')
    return start_point, last_vector_for_cu, last_vector_for_lyap

def find_min_ang(params, start_point, first_vector_for_cu, first_vector_for_lyap, count_of_trajectory, skip_var_num, points_main, points_skip, vectors_cu, vectors_ss):
    start = time.time()
    angles, min_ang = [None] * count_of_trajectory, 10
    der_norm = first_vector_for_lyap
    der_norm1 = first_vector_for_cu
    p1 = 0
    p2 = 0
    for i in range(count_of_trajectory):
        makeStepVar(der_norm, dimension, ShimizuX3_3D_flow_var, params, step, start_point)
        makeStepVar(der_norm1, dimension, ShimizuX3_3D_flow_var_trans_rev, params, step, start_point)
        makeStep(start_point, dimension, ShimizuX3_3D_flow, params, step)
        points_main[i] = start_point.copy()
        new_norm = np.linalg.norm(der_norm)
        new_norm1 = np.linalg.norm(der_norm1)
        der_norm = der_norm / new_norm
        der_norm1 = der_norm1 / new_norm1
        vectors_cu[i] = der_norm1.copy()
        p2 += m.log(new_norm)
        p1 += m.log(new_norm1)
    print('trans lyap', p1 / integrate_time)
    print('direct lyap', p2 / integrate_time)
    last_lyap_vector = der_norm
    last_point = start_point
    last_vector = der_norm1
    for i in range(skip_var_num):
        makeStep(start_point, dimension, ShimizuX3_3D_flow, params, step)
        points_skip[i] = start_point
    der = [eps, 0, 0]
    der_norm = der / np.linalg.norm(der)
    for i in range(skip_var_num):
        makeStepVar(der_norm, dimension, ShimizuX3_3D_flow_var_rev, params, step, points_skip[skip_var_num - 1 - i])
        new_norm = np.linalg.norm(der_norm)
        der_norm = der_norm / new_norm
    p3 = 0
    for i in range(count_of_trajectory):
        vectors_ss[count_of_trajectory - 1 - i] = der_norm.copy()
        angles[i] = abs(np.pi / 2 - np.arccos(
            np.dot(vectors_cu[count_of_trajectory - 1 - i], vectors_ss[count_of_trajectory - 1 - i])))
        makeStepVar(der_norm, dimension, ShimizuX3_3D_flow_var_rev, params, step, points_main[count_of_trajectory - 1 - i])
        p3 += m.log(np.linalg.norm(der_norm))
        der_norm = der_norm / np.linalg.norm(der_norm)
    print('back time lyap', p3 / integrate_time)
    # print(len(points_for_attr1[2]))
    # print(len(points_for_attr2[2]))
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(points_for_attr1[0], points_for_attr1[1], points_for_attr1[2], s=2, c='orange')
    # ax.scatter(points_for_attr2[0], points_for_attr2[1], points_for_attr2[2], s=5, c='blue')
    # plt.show()
    #vectors_ss.reverse()
    # for i in range(count_of_trajectory):
    #     angles[i] = abs(np.pi / 2 - np.arccos(np.dot(vectors_cu[i], vectors_ss[i])))
    for angle in angles:
        if min_ang > angle:
            min_ang = angle
    end = time.time()
    print(f'time for num: {end-start}')
    return angles, min_ang, last_point, last_vector, last_lyap_vector

def draw_continuity_cloud(axCU, axSS, skip, points, vectors_ss, vectors_cu, points_len):
    angles, dist = [], []
    # vectors_ss = np.flip(vectors_ss, 0)
    print(points_len/skip)
    for i in range(int(points_len / skip)):
        for j in range(int(points_len / skip)):
            if j <= i:
                continue
            angles.append(np.arccos(np.dot(vectors_ss[i*skip], vectors_ss[j*skip])/(np.linalg.norm(vectors_ss[i*skip])*np.linalg.norm(vectors_ss[j*skip]))))
            dist.append(m.dist(points[i*skip], points[j*skip]))
    print(len(dist), len(angles))
    axSS.scatter(dist, angles, c='#1F77B4', s=0.05)

    angles = []
    for i in range(int(points_len / skip)):
        for j in range(int(points_len / skip)):
            if j <= i:
                continue
            angles.append(np.arccos(np.dot(vectors_cu[i*skip], vectors_cu[j*skip])/(np.linalg.norm(vectors_cu[i*skip])*np.linalg.norm(vectors_cu[j*skip]))))
            # dist.append(m.dist(points[i*skip], points[j*skip]))
    print(len(dist), len(angles))
    axCU.scatter(dist, angles, c='#1F77B4', s=0.05)


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
fig.suptitle('CU', fontsize=20)
draw_continuity_cloud(axCu, axSs, skip_for_phase, points_main, vectors_ss, vectors_cu, int(integrate_time / step / num))
print(f'min angle: {min_ang}')

figPP = plt.figure()
axPP = figPP.add_subplot(projection='3d')
pT=points_main.T
axPP.plot(pT[0], pT[1], pT[2]) # отрисовка фазового портрета
plt.show()
# with open("file.txt", "w") as output:
#     output.write(str(all_angles))
