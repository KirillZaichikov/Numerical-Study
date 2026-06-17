from Step import *
import numpy as np
from model import *


# ПАРАМЕТРЫ ДЛЯ LERMANMAP S2
params = np.array([0.05, 0.7, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 ПРОСТО
# params = np.array([0.05, 0.041, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 ТУТ NAN
# params = np.array([0.0002, 0.5, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 ПАРА КРУГОВ
# params = np.array([0.0004, 0.5, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 ПАРА КРУГОВ
# params = np.array([0.00109, 0.85, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 ПАРА КРУГОВ
# params = np.array([0.00165, 0.53, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 Уже страшный хаос
# params = np.array([0.001, 0.055, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 Наш квазиаттрактор
# params = np.array([0.001, 0.02, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 Наш квазиаттрактор
# params = np.array([0.00001, 0.006, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 Уже страшный хаос
# params = np.array([0.0009, 0.4, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 РЕЗОНАНС КРУГОВ
# params = np.array([0.00056, 0.42, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 РЕЗОНАНС КРУГОВ
# params = np.array([0.00088, 0.66, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 КРАСИВЫЙ ХАОС

# ПАРАМЕТРЫ ДЛЯ LERMANMAP S3
# params = np.array([0.05, 0.2, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS2 ПРОСТО
# params = np.array([0.0002, 0.5, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  4 КРУГА
# params = np.array([0.00038, 0.52728, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  РЕЗОНАНС
# params = np.array([0.00045, 0.500, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  РЕЗОНАНС
# params = np.array([0.000484108, 0.39208, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  Мультистабильность
# params = np.array([0.000490921, 0.39208, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  Мультистабильность
# params = np.array([0.0002, 0.8, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  TestLyap
# params = np.array([0.0001, 0.8, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  TestLyap
# params = np.array([0.014, 0.048, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS3  TestLyap

# ПАРАМЕТРЫ ДЛЯ LERMANMAP S1
# params = np.array([0.0014, 0.195, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS1 
# params = np.array([0.006, 0.45, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS1 ПЕРВЫЙ ФАЗОВЫЙ
# params = np.array([0.03, 0.41, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS1 ВТОРОЙ ФАЗОВЫЙ
# params = np.array([0.15, 0.45, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMS1 ТРЕТИЙ ФАЗОВЫЙ
# params = np.array([0.02, 0.045, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRMTest
# params = np.array([-0.34, 0, 0], dtype=np.longdouble) # eps u v testModel
# params = np.array([1/3], dtype=np.longdouble) # eps alpha beta p smaleWilliams
# params = np.array([0.3, 0.4, 0.5, 5, 1], dtype=np.longdouble) # mu eps alpha beta p lerFRMwithMu_3D_map
# params = np.array([0.387, 0.56, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRM_3D_map

model_type = 'map'
diffFunc = lerFRM_3D_map
# diffFuncRev = lerFRM_3D_map_rev
diffFuncVarF = None # или None
dimension = 3

step = 1
skip_time = 10000
integrate_time = 3000000
skip_for_phase = 1
# initial_point = np.zeros(dimension)
# initial_point[0] = -0.7
# initial_point[1] = 0.1
# initial_point[2] = -0.7
initial_point = [[0.001,0.001,0.001]]
# initial_point = [[0.001,0.001,0.001],[0.5,0.1,0.0]]
# initial_point[0] = -params[0] + eps
# initial_point[1] = -params[1] / params[0]

# initial_point = np.array([0.001, 0.001, 0.001])
# initial_point = np.array([-0.0316478,  -0.06491352, -2.92705166]) # LermanMap Неподвижная точка (абсолютно неуст) [0.05, 0.43, 0.2, 1]
# initial_point = np.array([0.03395866, 0.01807846, 0.87946063]) # LermanMap точка периода 5 (седло 1 2) [0.05, 0.43, 0.2, 1]

crossection = 1 # for Lorenz3d 25 (z) for Lorenz4d 30/40 (y) for Shimizu3X 1.4(z)
eps = 1e-6 # for Lorenz3d 0.4 1 # for Shimizu3X 0.01

"""for lyapunov"""
lyap_num = 3
skip_var_time = 10000

"""for kneadings"""
seq_len = 500

"""for angles"""
with_clouds = False


if model_type == 'map':
    step = 1
    makeStep = MapStep
    makeStepVar = MapStepVarMat
