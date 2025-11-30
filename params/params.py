import numpy as np
from utils.integrator import *
from models.model import *


# params = np.array([1, -1, 1.8, 0.8, 0.5]) # for gonchenkoMap
# params = np.array([1.117, 0.5])
# params = np.array([10, 28, 8/3]) # for 3DLorenz
# params = np.array([1.367, 0.4, 0.2]) # for Shimizu3X
# params = np.array([10, 25, 8/3, 7]) # for 4dLorenz # 4.186139965, 44.08486,
# params = np.array([1, 5])
# params = np.array([-0.1, 5.1, 0])
# params = np.array([1.371672077922078, 0.43540322580645163, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 0
# params = np.array([1.3819673427222132, 0.4355835067637877, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 1
# params = np.array([0.68, 0.61, 0.2], dtype=np.longdouble) # alpha lamda B HIGH PRIORITY 2
params = np.array([4.3281, 39, 8/3, 7]) # alpha lamda
# params = np.array([1.0501456818181818, 0.6030381935483871, 0.2]) # alpha lamda B for Shimizu3X
params = np.array([0.05, 0.3, 0.2, 1], dtype=np.longdouble) # eps alpha beta p lerFRM
params = np.array([0.3, 0.4, 0.5, 5, 1], dtype=np.longdouble) # mu eps alpha beta p lerFRMwithMu_3D_map
model_type = 'map'
diffFunc = lerFRMwithMu_3D_map
diffFuncVarF = None # или None
dimension = 3

step = 0.01
skip_time = 200000
integrate_time = 100000
skip_for_phase = 1
initial_point = np.zeros(dimension)
# initial_point[0] = -params[0] + eps
# initial_point[1] = -params[1] / params[0]
initial_point[0] = 0.001
initial_point[1] = 0.001
initial_point[2] = 0.0001

crossection = 1 # for Lorenz3d 25 (z) for Lorenz4d 30/40 (y) for Shimizu3X 1.4(z)
eps = 0.01 # for Lorenz3d 0.4 1 # for Shimizu3X 0.01

"""for lyapunov"""
lyap_num = 3
skip_var_time = 200

"""for kneadings"""
seq_len = 500

makeStep = dverkStep
makeStepVar = dverkStepVarMat
# diffFuncVar = eval(f"{diffFunc}" + "var")
if model_type == 'map':
    step = 1
    makeStep = MapStep
    makeStepVar = MapStepVarMat
