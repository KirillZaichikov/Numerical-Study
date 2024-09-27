import time
import math as m
from utils.linal import *
from params.params import *
from utils.integrator import *

def main():
    der = np.identity(dimension)
    new_norm = np.zeros(lyap_num)
    P = [0 for i in range(dimension)]

    start = time.time()
    for i in range(int(skip_time/step)):
        makeStep(initial_point, dimension, diffFunc, params, step)
    print("***********TIME for SKIP LYAPUNOV:", time.time() - start)

    start = time.time()
    for i in range(int(skip_var_time/step)):
        for j in range(lyap_num):
            makeStepVar(der[j], dimension, diffFuncVarF, params, step, initial_point)
            new_norm[j] = np.linalg.norm(der[j])
            der[j] = der[j] / new_norm[j]
        ortVecs(der, dimension, lyap_num)
        makeStep(initial_point, dimension, diffFunc, params, step)


    for i in range(int(integrate_time/step)):
        for j in range(lyap_num):
            makeStepVar(der[j], dimension, diffFuncVarF, params, step, initial_point)
            new_norm[j] = np.linalg.norm(der[j])
            der[j] = der[j] / new_norm[j]
            P[j] = P[j] + m.log(new_norm[j])
        ortVecs(der, dimension, lyap_num)
        makeStep(initial_point, dimension, diffFunc, params, step)
    print("***********TIME for SKIP_VAR + MainTraj:", time.time() - start)

    for i in range(lyap_num):
        print(f"L{i+1}: ", P[i] / integrate_time)

if __name__ == "__main__":
    main()