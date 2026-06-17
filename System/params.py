import numpy as np
from integrator import *
from system import *

# params = np.array([1, 1, 0.3, 0.25, 0.0, 1, 0.8]) # ro, m, i1, i2, a0, g, en
params = np.array([1, 1, 0.3, 0.25, 0.55, 1, 5]) # ro, m, i1, i2, a0, g, en
# Консервативный если энергия больше 1.55
# Или чтобы локально и консервативно, то инерции одинаковые, а en меньше 1

model_type = 'flow'
diffFunc = omegaGamma_sys
# diffFunc = main_sys
diffPoincare = main_sys_poincare # если сечение по phi
# diffPoincare = omegaGamma_sys_poincare # если сечение по omega1
diffFuncRev = None
diffFuncVarF = None # или None
dimension = 5
dimension_after_replace = 3

step = 0.00001
skip_time = 0
integrate_time = 100
skip_for_phase = 1

initial_point = np.zeros(dimension_after_replace)
# initial_point[0] = 0.5
# initial_point[1] = 3*np.pi/2
# initial_point[2] = (np.pi/4)-0.01

initial_point[0] = 0.01
initial_point[1] = 0.5
initial_point[2] = 0.01

crossection = 0.5
# crossection = 0
DISCR = 100
iter_num = 5

eps = 1e-2

"""for lyapunov"""
lyap_num = 3
skip_var_time = 2000

"""for kneadings"""
seq_len = 500

"""for angles"""
with_clouds = False

makeStep = dverkStep
makeStepVar = dverkStepVarMat
