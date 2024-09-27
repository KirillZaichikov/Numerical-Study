import numpy as np
from utils.integrator import *
from models.model import *

model_type = 'flow'
# params = np.array([1, -1, 1.8, 0.8, 0.5]) # for gonchenkoMap
# params = np.array([1.117, 0.5])
# params = np.array([10, 28, 8/3]) # for 3DLorenz
# params = np.array([4.74, 34.0751, 8/3, 7]) # for 4dLorenz
params = np.array([1.367, 0.4, 0.2]) # for Shimizu3X
diffFunc = ShimizuX3_3D_flow
diffFuncVarF = None #f"{diffFunc}_var"
dimension = 3
delta_params = 1 # for Lorenz3d 25 (z) for Lorenz4d 30/40 (y) for Shimizu3X 1.4(z)
eps = 0.01 # for Lorenz3d 0.4 1 # for Shimizu3X 0.01
step = 0.001
skip_time = 0
integrate_time = 1000
skip_for_phase = 10
initial_point = np.zeros(dimension)
initial_point[0] = 1e-40

"""for lyapunov"""
lyap_num = 3
skip_var_time = 100

makeStep = dverkStep
makeStepVar = dverkStepVarMat
# diffFuncVar = eval(f"{diffFunc}" + "var")
if model_type == 'map':
    step = 1
    makeStep = MapStep
    makeStepVar = MapStepVarMat