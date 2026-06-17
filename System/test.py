import numpy as np
from params import *

def norm_solution(start_point):
    if start_point[2] > 2*np.pi:
        start_point[2] -= 2*np.pi
    elif start_point[2] < 0:
        start_point[2] += 2*np.pi

def mg_to_oop(start_point):
    '''
    MG -> omega1, omega2, phi
    '''
    omega1, omega2, _g1, _g2, _g3 = start_point
    _phi = np.atan2(_g1, _g2)
    print(_phi)
    point = np.array([omega1, omega2, _phi])
    norm_solution(point)
    # ТУТ ОШИБОК НЕ ВОЗНИКАЕТ
    if point[2] > 2*np.pi or point[2] < 0:
        print('ERROR(mg_to_oop)', point)
        print(np.sqrt(-1))
    return point.copy()

def mg_to_dpt(start_point, params):
    '''
    MG -> dzeta, phi, teta
    '''
    omega1, omega2, _g1, _g2, _g3 = start_point
    _phi = np.atan2(_g1, _g2)
    _teta = np.acos(_g3)
    _ksi = np.atan2(omega2, omega1)
    point = np.array([_ksi+_phi,_phi,_teta])
    norm_solution(point)
    # ТУТ ОШИБОК НЕ ВОЗНИКАЕТ (КАЖЕТСЯ ВОЗНИКАЕТ НО НЕ ОБРАБАТЫВАЕТСЯ)
    if _phi>2*np.pi or _phi<0 or point[2] < 0 or point[2] > 2*np.pi:
        print('ERROR(mg_to_dpt)', point, _phi)
    return point.copy()

params = np.array([1, 1, 0.3, 0.25, 0.55, 1, 5]) # ro, m, i1, i2, a0, g, en
gamma1 = -0.2
gamma2 =  0.3
gamma3 =  0.9
norm = np.linalg.norm([gamma1,gamma2,gamma3])
print(mg_to_oop([1,2, gamma1/norm, gamma2/norm, gamma3/norm]))
print(mg_to_dpt([1,2, gamma1/norm, gamma2/norm, gamma3/norm], params))