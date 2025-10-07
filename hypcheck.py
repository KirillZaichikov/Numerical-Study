from utils.integrator import *
from utils.linal import *
import math as m
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

dimension = 3
step = 0.001
eps = 1

skip_time = 0
skip_var_time = 60
integrate_time = 250
num = 1

initial_point = np.array([1e-50, 0, 0])
params = np.array([1.18, 0.4])

def Shimizu_3D_flow(state, res, params, H) -> None:
    res[0] = state[1]
    res[1] = - params[0] * state[1] - state[0] * state[2] + state[0]
    res[2] = - params[1] * state[2] + pow( state[0] , 2.0 )

def Shimizu_3D_flow_var(state, res, params, stateOld):
   res[0] =                                                 state[1] *  1 
   res[1] = state[0] * ( 1 - stateOld[2] )                + state[1] * ( - params[0] ) + state[2] * ( - stateOld[0] ) 
   res[2] = state[0] * ( 2.0 * pow( stateOld[0] , 1.0 ) )                              + state[2] * ( - params[1] ) 

def Shimizu_3D_flow_var_trans(state, res, params, stateOld):
   res[0] =                 state[1] * ( 1 - stateOld[2] ) + state[2] * ( 2.0 * pow( stateOld[0] , 1.0 ) )     
   res[1] = state[0] *  1 + state[1] * ( - params[0] )
   res[2] =                 state[1] * ( - stateOld[0] )   + state[2] * ( - params[1] ) 

def Shimizu_3D_flow_var_rev(state, res, params, stateOld):
   res[0] = -1 * ( state[0] * ( 0 ) +state[1] * ( 1 ) +state[2] * ( 0 )  )
   res[1] = -1 * ( state[0] * ( 1 - stateOld[2] ) +state[1] * ( - params[0] ) +state[2] * ( - stateOld[0] )  )
   res[2] = -1 * ( state[0] * ( 2.0 * pow( stateOld[0] , 1.0 ) ) +state[1] * ( 0 ) +state[2] * ( - params[1] )  )

def Shimizu_3D_flow_var_trans_rev(state, res, params, stateOld):
   res[0] = -1 * (state[1] * ( 1 - stateOld[2] ) + state[2] * ( 2.0 * pow( stateOld[0] , 1.0 ) ) )
   res[1] = -1 * (state[0] *  1  + state[1] * ( - params[0] ) )
   res[2] = -1 * (state[1] * ( - stateOld[0] ) + state[2] * ( - params[1] ) )

def func_for_get_to_attractor(params, start_point, skip_num, skip_var_num):
    der = [eps, 0, 0]
    der_norm = der / np.linalg.norm(der)
    for i in range(skip_num):
        dverkStep(start_point, dimension, Shimizu_3D_flow, params, step)
    for i in range(skip_var_num):
        dverkStepVarMat(der_norm, dimension, Shimizu_3D_flow_var, params, step, start_point)
        dverkStep(start_point, dimension, Shimizu_3D_flow, params, step)
        new_norm = np.linalg.norm(der_norm)
        der_norm = der_norm / new_norm

    last_vector_for_lyap =  der_norm
    return start_point, last_vector_for_lyap

def find_min_ang(params, start_point, first_vector_for_lyap, count_of_trajectory, skip_var_num, points_main, points_skip, vectors_u, vectors_cs):
    angles, min_ang = [None] * count_of_trajectory, 10
    der_norm = first_vector_for_lyap
    p1 = 0
    for i in range(count_of_trajectory):
        dverkStepVarMat(der_norm, dimension, Shimizu_3D_flow_var, params, step, start_point)
        dverkStep(start_point, dimension, Shimizu_3D_flow, params, step)
        points_main[i] = start_point.copy()
        new_norm = np.linalg.norm(der_norm)
        der_norm = der_norm / new_norm
        vectors_u[i] = der_norm.copy()
        p1 += m.log(new_norm)
    print('direct lyap', p1 / integrate_time)
    last_lyap_vector = der_norm
    last_point = start_point
    for i in range(skip_var_num):
        dverkStep(start_point, dimension, Shimizu_3D_flow, params, step)
        points_skip[i] = start_point
    der = [eps, 0, 0]
    der_norm = der / np.linalg.norm(der)
    for i in range(skip_var_num):
        dverkStepVarMat(der_norm, dimension, Shimizu_3D_flow_var_trans, params, step, points_skip[skip_var_num - 1 - i])
        new_norm = np.linalg.norm(der_norm)
        der_norm = der_norm / new_norm
    p3 = 0
    for i in range(count_of_trajectory):
        vectors_cs[count_of_trajectory - 1 - i] = der_norm.copy()
        angles[i] = abs(np.pi / 2 - np.arccos(np.dot(vectors_u[count_of_trajectory - 1 - i],
                                                     vectors_cs[count_of_trajectory - 1 - i])))
        dverkStepVarMat(der_norm, dimension, Shimizu_3D_flow_var_trans, params, step, points_main[count_of_trajectory - 1 - i])
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
    return angles, min_ang, last_point, last_lyap_vector


min_ang = 10
all_angles = []
last_point, last_vector_for_lyap = func_for_get_to_attractor(params, initial_point, int(skip_time / step), int(skip_var_time / step))
points_main = np.zeros(shape=(int(integrate_time / step / num), dimension), dtype=np.float32)
points_skip = np.zeros(shape=(int(skip_var_time / step), dimension), dtype=np.float32)
vectors_cu, vectors_ss = np.zeros(shape=(int(integrate_time / step / num), dimension), dtype=np.float32), np.zeros(shape=(int(integrate_time / step / num), dimension), dtype=np.float32)

for i in range(int(num)):
    angles, angle, last_point, last_vector_for_lyap = find_min_ang(params, last_point, last_vector_for_lyap,
                                                        int(integrate_time / step / num), int(skip_var_time/step), points_main, points_skip, vectors_cu, vectors_ss)
    if min_ang > angle:
        min_ang = angle
    all_angles.append(angles)

print(f'min angle: {min_ang}')

angles_np = np.array(angles)
print(angles_np)
print(f'max angle: {angles_np.max()}')
figPP = plt.figure()
axPP = figPP.add_subplot(projection='3d')
pT=points_main.T
# axPP.plot(pT[0], pT[1], pT[2]) # отрисовка фазового портрета

# threshold = 0.01
# norm = Normalize(vmin=angles_np.min(), vmax=threshold)
# scatter = axPP.scatter(pT[0][::-1], pT[1][::-1], pT[2][::-1], c=angles_np, cmap='viridis', s=0.01, norm=norm)
# cbar = plt.colorbar(scatter)
# cbar.set_label('Angle values')

pos = []
neg = []
for i, a in enumerate(angles_np):
    if a < 0.001:
        neg.append(points_main[i])
    else : 
        pos.append(points_main[i])
posN = np.array(pos)
negN = np.array(neg)
posN = posN.T
negN = negN.T
# if len(posN) != 0 :
#     axPP.scatter(posN[0], posN[1], posN[2], c="orange", s=0.005)
if len(negN) != 0 :
    axPP.scatter(negN[0], negN[1], negN[2], c="b", s=0.02)

plt.show()
