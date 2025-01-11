# import numpy as np
# import math as m
# import matplotlib.pyplot as plt
import time
from params.params import *
from utils.integrator import *
from models.model import *


# for phase portret
def main():
    print("***********TIME for integrate: ",integrate_time)

    start = time.time()
    eps = 1e-8
    last_was_ret = False
    while eps > 1e-20:
        find_hom = False
        initial_point = np.array([1e-20, 0, 0, 0])
        for i in range(int(integrate_time/step)):
            makeStep(initial_point, dimension, diffFunc, params, step)
            if initial_point[0] < 0:
                print("here hom", initial_point[0])
                find_hom = True
                break
        if not find_hom:
            params[0] = params[0] + eps
            eps = eps / 2
            params[0] = params[0] - eps
            last_was_ret = True
        elif last_was_ret:
            eps = eps / 2
            params[0] = params[0] - eps
            last_was_ret = False
        else:
            params[0] = params[0] - eps
            last_was_ret = False
        print(params[0])
    end = time.time()
    print("***********TIME:", end-start)


if __name__ == "__main__":
    main()