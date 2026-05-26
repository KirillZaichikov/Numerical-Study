'''
Тут пытаюсь считать углы, пока неудачно
'''

import numpy as np
import matplotlib.pyplot as plt
import math as m
from tqdm import tqdm
import time as t
import sys
sys.path.append(r'C:\Users\Kirill\Desktop\repos\Numerical-Study')
from utils.linal import *
from utils.integrator import *
from calcKuzFunctions import *


if __name__ == "__main__":
    x1, x2 = [], []

    # Начальная точка из хаоса
    M = np.array([176.39671645238903, -168.85257112196294, -94.877918535603939])
    gamma = np.array([0.69524430378189617, 0.40429993079272591, -0.59428864887075838])
    gamma = -gamma
    main_traj = np.array([*M, *gamma])

    # Параметры для вычислений
    step = 0.001
    time_skip = 0
    time_skip_clv = 30
    time = 100
    lyap_num = 6
    dimension = 6
    eps = 0.001

    # Параметры из хаоса
    params = [0.485, 2, 6, 7, 9, 4, 1, 752, 100] # d I1 I2 I3 a1 a2 h E g0
    # params = [0.2, 5, 6, 7, 9, 4, 1, 555, 100]
    d = params[0]
    I1, I2, I3 = params[1], params[2], params[3]
    a1, a2 = params[4], params[5]
    h, E, g0 = params[6], params[7], params[8]

    # Замена начальных условий из хаоса на значения кузнецова
    QC = np.array([[np.cos(d), np.sin(d), 0],
                 [-np.sin(d),np.cos(d),0],
                 [0,0,1]])
    M = QC @ M
    gamma = QC @ gamma

    # проверим интегралы
    En, gi = calc_integrals(M, gamma, d, I1, I2, I3, a1, a2, h, E, g0)

    iterNumCLP = int(time_skip_clv/step)
    iterNumMain = int(time/step)
    new_norm = np.zeros(dimension)
    P = [0 for i in range(dimension)]
    Pb = [0 for i in range(dimension)]
    slaveTrajectory = np.identity(dimension) * eps
    points_main = np.zeros(shape=(int(time / step), dimension), dtype=np.float32)
    points_skip = np.zeros(shape=(int(time_skip_clv / step), dimension), dtype=np.float32)
    vectors_ncu = np.zeros(shape=(int(time / step), dimension), dtype=np.float32)

    test1 = []
    test2 = []
    test2x = []
    main_traj = np.array([*M, *gamma])

    # Выход на аттрактор
    for i in tqdm(range (int(time_skip/step))):
        dverkStep(main_traj, dimension, sys, params, step)

    # Прогрев векторов
    for i in tqdm(range (int(iterNumCLP))):
        for j in range(lyap_num):
            slaveTrajectory[j] += main_traj

        dverkStep(main_traj, dimension, sys, params, step)

        for j in range(lyap_num):
            dverkStep(slaveTrajectory[j], dimension, sys, params, step)
            slaveTrajectory[j] -= main_traj
            # print(slaveTrajectory[j])

        ortVecs(slaveTrajectory, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(slaveTrajectory[j])
            slaveTrajectory[j] = (slaveTrajectory[j] / new_norm[j]) * eps

    # Считаем показатели в прямом времени и запоминаем точки и вектора
    for i in tqdm(range(iterNumMain)):
        for j in range(lyap_num):
            slaveTrajectory[j] += main_traj
        # test1.append(main_traj[0])
        # test2x.append(i)
        # test2.append(slaveTrajectory[0][0])

        dverkStep(main_traj, dimension, sys, params, step)

        points_main[i] = main_traj.copy()

        for j in range(lyap_num):
            dverkStep(slaveTrajectory[j], dimension, sys, params, step)
            slaveTrajectory[j] -= main_traj
            # if j==0:
            #     test2.append(slaveTrajectory[0][0])
            #     test2x.append(i+1)
        # print("dire", slaveTrajectory[0])
        ortVecs(slaveTrajectory, dimension, lyap_num)
        ### Раньше было в цикле возможно перестановка ни на что не влияет
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(slaveTrajectory[j])
            P[j] += m.log(new_norm[j] / eps)
            slaveTrajectory[j] = (slaveTrajectory[j] / new_norm[j]) * eps
        ###
        vectors_ncu[i] = slaveTrajectory[-1].copy()
    # print(vectors_ncu)

    # Запоминаем точки для прогрева в обратном времени
    for i in tqdm(range(iterNumCLP)):
        dverkStep(main_traj, dimension, sys, params, step)
        points_skip[i] = main_traj.copy()

    # slaveBackTime = np.array([eps, 0, 0, 0, 0, 0])
    slaveBackTime = np.identity(dimension) * eps
    for j in range(lyap_num):
        slaveBackTime[j] += points_skip[-1]
    # Прогрев в обратном времени
    for i in tqdm(range(2, iterNumCLP+1)): # Думаю здесь не надо делать последнюю итерацию (Вшил защиту в цикл)

        # if i == (iterNumCLP - 1): # Думаю так надо
        # for j in range(lyap_num):
        #     dverkStep(slaveBackTime[j], dimension, sysRev, params, step)
        # ortVecs(slaveBackTime, dimension, lyap_num)
        #     # break

        for j in range(lyap_num):
            dverkStep(slaveBackTime[j], dimension, sysRev, params, step)
            slaveBackTime[j] -= points_skip[-i]
        ortVecs(slaveBackTime, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(slaveBackTime[j])
            slaveBackTime[j] = (slaveBackTime[j] / new_norm[j]) * eps
        for j in range(lyap_num):
            slaveBackTime[j] += points_skip[-i]

    print(slaveBackTime)
    for j in range(lyap_num):
        dverkStep(slaveBackTime[j], dimension, sysRev, params, step)
    print(slaveBackTime)
    for j in range(lyap_num):
        slaveBackTime[j] -= points_main[-1]
    print(slaveBackTime)
    ortVecs(slaveBackTime, dimension, lyap_num)
    for j in range(lyap_num):
        new_norm[j] = np.linalg.norm(slaveBackTime[j])
        slaveBackTime[j] = (slaveBackTime[j] / new_norm[j]) * eps

    # Зачем то доделываю с последней итерацией
    # for j in range(lyap_num):
    #     slaveBackTime[j] -= points_main[-1]
    #     new_norm[j] = np.linalg.norm(slaveBackTime[j])
    #     slaveBackTime[j] = (slaveBackTime[j] / new_norm[j]) * eps
        # print(slaveBackTime[j])
    # ortVecs(slaveBackTime, dimension, lyap_num)

    test1b = []
    test2b = []
    test2xb = []

    min_angle = 10
    angles = []
    # print(slaveBackTime)
    # for i in range(iterNumMain-1):
    for i in tqdm(range(iterNumMain-1)): # Здесь думаю точно последнюю надо пропустить, нет вектора из которого вычесть
        # angle = abs(np.pi / 2 - np.arccos(np.dot(slaveBackTime[0]/ eps, vectors_ncu[iterNumMain - 1 - i]/eps)))
        # angles.append(angle)
        # if angle < min_angle:
        #     min_angle = angle
        # test1b.append(points_main[iterNumMain - 1 - i][0])
        for j in range(lyap_num):
            slaveBackTime[j] += points_main[iterNumMain - 1 - i]
            # if j==0:
            #     test2xb.append(iterNumMain - i)
            #     test2b.append(slaveBackTime[j][0])
            dverkStep(slaveBackTime[j], dimension, sysRev, params, step)
            # if j==0:
            #     test2xb.append(iterNumMain - 1 - i)
            #     test2b.append(slaveBackTime[j][0])
            slaveBackTime[j] -= points_main[iterNumMain - 2 - i]

        # print("back", slaveBackTime[0])
        ortVecs(slaveBackTime, dimension, lyap_num)
        for j in range(lyap_num):
            new_norm[j] = np.linalg.norm(slaveBackTime[j])
            slaveBackTime[j] = (slaveBackTime[j] / new_norm[j]) * eps
            Pb[j] += m.log(new_norm[j] / eps)


    for i in range(lyap_num):
        print(f"L{i+1}: ", P[i] / time)
    for i in range(lyap_num):
        print(f"Lb{i+1}: ", Pb[i] / time)


    # проверим интегралы
    En, gi = calc_integrals(M, gamma, d, I1, I2, I3, a1, a2, h, E, g0)

    # print('min angle', min_angle)
    # x = np.linspace(0, len(test1)-1, len(test1))
    # xb = np.linspace(2, len(test1b)+1, len(test1b))
    # fig, axes = plt.subplots()
    # axes.plot(x,test1)
    # axes.plot(test2x,test2)
    # test1b.reverse()
    # axes.plot(xb,test1b)
    # axes.plot(test2xb,test2b)

    plt.show()