import numpy as np

mas_for_poincare = np.zeros((10, 2))
mas_for_poincare[0] = [1, 1]
mas_for_poincare[4] = [1, 0]
mas_for_poincare[8] = [2, 2]
mas_for_poincare[9] = [3, 3]
# print(mas_for_poincare)
# print(mas_for_poincare[mas_for_poincare != 0])

mas_for_poincare = np.array([ [None] * 2 for i in range(10)])
mas_for_poincare[0] = [1, 1]
mas_for_poincare[4] = [1, 0]
mas_for_poincare[8] = [2, 2]
mas_for_poincare[9] = [3, 3]
print(mas_for_poincare)
# mas_for_poincare = mas_for_poincare[mas_for_poincare != None]
print(mas_for_poincare)
# mas_for_poincare = mas_for_poincare.reshape()
# print(mas_for_poincare)
# mas_for_poincare = mas_for_poincare.T
# print(mas_for_poincare[0])