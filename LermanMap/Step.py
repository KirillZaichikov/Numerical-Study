import numpy as np

def MapStep(val, dimension, Func, params, step) -> None:
    res = np.zeros(dimension)
    Func(val, res, params)
    for j in range(dimension):
        val[j] = res[j]
    # Func(val, res, params)
    # for j in range(dimension):
    #     val[j] = res[j]


def MapStepVarMat(val, dimension, FuncVar, params, step, mainTraj) -> None:
    res = np.zeros(dimension)
    FuncVar(val, res, params, mainTraj)
    for j in range(dimension):
        val[j] = res[j]