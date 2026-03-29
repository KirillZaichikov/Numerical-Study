import math as m
import numpy as np
import matplotlib.pyplot as plt

alpha = 1/3
params = [alpha]

def Rev(state,res,params):
    phi_, X_, Y_ = state
    alpha = params[0]
    
    if phi_ > 3*np.pi/2:
        phi = phi_ / 2
    else:
        phi = (phi_ + 2 * np.pi) / 2
    # if (phi > 2*m.pi):
    #     phi = phi - 2 * m.pi
    # elif (phi < 0):
    #     phi += 2 * m.pi
    X = (X_ - (1/2)*np.cos(phi)) / alpha
    Y = (Y_ - (1/2)*np.sin(phi)) / alpha

    res[0] = phi
    res[1] = X
    res[2] = Y

def Dir(state,res,params):
    phi, X, Y = state
    alpha = params[0]

    phi_ = 2 * phi
    if (phi_ > 2*m.pi):
        phi_ = phi_ - 2 * m.pi
    elif (phi_ < 0):
        phi_ += 2 * m.pi
    X_ = alpha * X + (1/2)*np.cos(phi)
    Y_ = alpha * Y + (1/2)*np.sin(phi)

    res[0] = phi_
    res[1] = X_
    res[2] = Y_

iter_num=10
y1, y2 = [], []
initial_point = np.array([5.250885, 0.27248168, -0.27191988])
res_point = np.zeros(3, dtype=float)
for i in range(iter_num):
    Dir(initial_point, res_point, params)
    print(res_point)
    # for j in range(3):
    #     initial_point[j] = res_point[j]
    initial_point = res_point.copy()
    y1.append(res_point.copy())
print("rev start")
for i in range(iter_num):
    y2.append(res_point.copy())
    initial_point = res_point.copy()
    print(res_point)
    Rev(initial_point, res_point, params)

t = np.arange(0, iter_num)
yn1 = np.array(y1)
yn2 = np.array(y2)
y1 = yn1.T
y2 = yn2.T

fig, ax = plt.subplots()
ax.plot(t, y1[1])
ax.plot(t, y2[1][::-1])
fig1, ax1 = plt.subplots()
ax1.plot(y1[1],y1[2],linestyle="", marker = '.', ms=5)
ax1.plot(y2[1],y2[2],linestyle="", marker = '.', ms=5)
plt.show()