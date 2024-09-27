import numpy as np

def dotProduct(firstVec, secondVec, length):
    res = 0
    for i in range(length):
        res += firstVec[i] * secondVec[i]
    return res


def ortVecs(vecs, dimension, numVecs):
    projCoeff = 0
    projSum = np.zeros(dimension)

    for i in range(1,numVecs):
        for j in range(i):
            projCoeff = dotProduct(vecs[j], vecs[i], dimension) / dotProduct(vecs[j], vecs[j], dimension)
            for k in range(dimension):
                projSum[k] += projCoeff * vecs[j][k]
        
        for j in range(dimension):
            vecs[i][j] -= projSum[j]
            projSum[j] = 0
