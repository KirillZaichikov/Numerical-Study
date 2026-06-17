import time
import math as m
from utils.linal import *
from params.params import *
from utils.integrator import *
# from numba import jit

# @jit(nopython=True, cache=True)
def main(initial_point):
    der = np.identity(dimension) * eps
    # print(der)
    new_norm = np.zeros(lyap_num)
    P = [0 for i in range(dimension)]

    # start = time.time()
    print("SKIP_TIME, STEP", skip_time, step)
    for i in range(int(skip_time/step)):
        makeStep(initial_point, dimension, diffFunc, params, step)
    # print("***********TIME for SKIP LYAPUNOV:", time.time() - start)

    # start = time.time()
    print("SKIP_VAR, STEP", skip_var_time, step)
    for i in range(int(skip_var_time/step)):
        for j in range(lyap_num):
            der[j] += initial_point
            makeStep(der[j], dimension, diffFunc, params, step)
            # print(der[0])
        makeStep(initial_point, dimension, diffFunc, params, step)
        # print(initial_point)
        for j in range(lyap_num):
            # if np.linalg.norm(der[j]) >= np.pi: print("ERROR_WARM")
            der[j] -= initial_point
        ortVecs(der, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(der[j])
            der[j] = (der[j] / new_norm[j])*eps

    for i in range(int(integrate_time/step)):
        for j in range(lyap_num):
            der[j] += initial_point
            makeStep(der[j], dimension, diffFunc, params, step)
        makeStep(initial_point, dimension, diffFunc, params, step)
        for j in range(lyap_num):
            # if new_norm[j] >= np.pi: print("ERROR_MAIN")
            der[j] -= initial_point
        ortVecs(der, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(der[j])
            der[j] = (der[j] / new_norm[j])*eps
            P[j] += m.log(new_norm[j]/eps)
        
    # print("***********TIME for SKIP_VAR + MainTraj:", time.time() - start)

    for i in range(lyap_num):
        print(f"L{i+1}: ", P[i] / integrate_time)


if __name__ == "__main__":
    start_point = np.array(initial_point)
    for i in start_point:
        print(i)
        main(i)