import numpy as np
import matplotlib.pyplot as plt

a = []
with open("matrix.txt","r") as file:
    for line in file.readlines():
        tmp = line[:-2].split(' ')
        for i in range(len(tmp)):
            tmp[i] = float(tmp[i])
        a.append(tmp)

mas = np.array(a)
coord = mas.T
# print(coord[0])

fig_phase = plt.figure()
ax_phase = fig_phase.add_subplot(projection='3d')
ax_phase.scatter(coord[0], coord[1], coord[3], s=0.1, c='black')
plt.show()
