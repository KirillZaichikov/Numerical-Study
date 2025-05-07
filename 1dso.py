import matplotlib.pyplot as plt
import numpy as np

mu_ = np.linspace(1.53577, 1.59377, 1000)#0.104
A_ = np.linspace(-1.152, -1.0511, 1000)
C_ = 1.3
nu_ = 0.8#0.964
skip_time = 0
time = 100
times = 1
skip_time = 0
x_ = 0.0
lx, ly = [], []

def map(x, A, C, nu, mu):
    return abs(-mu+A*x**nu+C*x**(2*nu))

def map_der(x, A, C, nu, mu):
    if -mu+A*x**nu+C*x**(2*nu) > 0:
        return A*nu * x **(nu-1) + 2 * nu * C * x ** (2*nu-1)
    else:
        return -(A*nu * x **(nu-1) + 2 * nu * C * x ** (2*nu-1))
    
points_pos = []
points_neg = []
points_nan = []
points_cvasi = []
for A in A_:
    print('one line is ready')
    for mu in mu_:
        # print(A, mu)
        x_ = 0
        min_deriv = 100
        cvasi = False
        for j in range(skip_time):
            x_ = map(x_, A, C_, nu_, mu)
        for j in range(time):
            x_ = map(x_, A, C_, nu_, mu)
            deriv = map_der(x_, A, C_, nu_, mu)
            
            # if abs(deriv) < min_deriv and time > 100:
            #     cvasi = True
            #     min_deriv = abs(deriv)
            if abs(deriv) < min_deriv:
                min_deriv = abs(deriv)
        # if cvasi:
        #     points_cvasi.append([A, mu])
        if min_deriv > 0.002:
            points_pos.append([A, mu])
        elif min_deriv == 'inf':
            points_nan.append([A, mu])
        else:
            # print(deriv, A, mu)
            points_neg.append([A, mu])

pos = np.array(points_pos).T
neg = np.array(points_neg).T
nan = np.array(points_nan).T
# cva = np.array(points_cvasi).T

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(mu_[0],mu_[-1])
ax.set_ylim(A_[0],A_[-1])

# ax.scatter(pos[1], pos[0], s=0.1, c='white')
# ax.scatter(cva[1], cva[0], s=0.05, c='blue')
ax.scatter(neg[1], neg[0], s=0.1, c='black')
if len(points_nan) != 0:
    ax.scatter(nan[1], nan[0], s=0.1, c='gray')
plt.show()