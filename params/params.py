import numpy as np
from utils.integrator import *
from models.model import *

model_type = 'flow'
# params = np.array([1, -1, 1.8, 0.8, 0.5]) # for gonchenkoMap
# params = np.array([1.117, 0.5])
# params = np.array([10, 28, 8/3]) # for 3DLorenz
# params = np.array([1.367, 0.4, 0.2]) # for Shimizu3X
params = np.array([10, 25, 8/3, 7]) # for 4dLorenz # 4.186139965, 44.08486,
# params = np.array([1, 5])
# params = np.array([-0.1, 5.1, 0])
diffFunc = Lorenz_4D_flow
diffFuncVarF = Lorenz_4D_flow_var # или None
dimension = 4
delta_params = 40 # for Lorenz3d 25 (z) for Lorenz4d 30/40 (y) for Shimizu3X 1.4(z)
eps = 0.01 # for Lorenz3d 0.4 1 # for Shimizu3X 0.01
step = 0.001
skip_time = 0
integrate_time = 100
skip_for_phase = 1
initial_point = np.zeros(dimension)
# initial_point[0] = -params[0] + eps
# initial_point[1] = -params[1] / params[0]
initial_point[0] = 0.001
initial_point[1] = 0
initial_point[2] = 0

"""for lyapunov"""
lyap_num = 3
skip_var_time = 50

"""for kneadings"""
seq_len = 500

makeStep = dverkStep
makeStepVar = dverkStepVarMat
# diffFuncVar = eval(f"{diffFunc}" + "var")
if model_type == 'map':
    step = 1
    makeStep = MapStep
    makeStepVar = MapStepVarMat