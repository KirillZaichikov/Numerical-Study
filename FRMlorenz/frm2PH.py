import matplotlib.pyplot as plt
import numpy as np
import math as m

LYAP_TRES = 0.001
DERI_TRES = 0.005

mu_ = np.linspace(-0.3, 3, 500)#0.104
A_ = np.linspace(-1.5, 1.1, 500)
C_ = 1.3
nu_ = 1#0.964
time = 1000
times = 1
skip_time = 200
x_ = 0.0
lx, ly = [], []

def map(x, A, C, nu, mu):
    return abs(-mu+A*x**nu+C*x**(2*nu))

def map_der(x, A, C, nu, mu):
    if -mu+A*x**nu+C*x**(2*nu) > 0:
        return A*nu * x **(nu-1) + 2 * nu * C * x ** (2*nu-1)
    else:
        return -(A*nu * x **(nu-1) + 2 * nu * C * x ** (2*nu-1))
    
points_ph = []
points_neg = []
points_nan = []
points_cvasi = []
for A in A_:
    print('one line is ready')
    for mu in mu_:
        # print(A, mu)
        lyap = 0
        x_ = 0
        min_deriv = 100
        cvasi = False
        for j in range(skip_time):
            x_ = map(x_, A, C_, nu_, mu)
        for j in range(time):
            x_ = map(x_, A, C_, nu_, mu)
            deriv = map_der(x_, A, C_, nu_, mu)
            try:
                if type(lyap) is not str:
                    lyap += m.log(deriv)
            except ValueError:
                lyap = 'inf'
            
            if abs(deriv) < min_deriv:
                min_deriv = abs(deriv)

        if min_deriv == 'inf' or lyap=='inf':
            points_nan.append([A, mu])
        elif min_deriv > DERI_TRES and lyap>LYAP_TRES:
            points_ph.append([A, mu])
        elif min_deriv < DERI_TRES and lyap>LYAP_TRES:
            points_cvasi.append([A, mu])
        elif lyap < LYAP_TRES:
            points_neg.append([A, mu])

ph = np.array(points_ph).T
neg = np.array(points_neg).T
nan = np.array(points_nan).T
cva = np.array(points_cvasi).T

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(mu_[0],mu_[-1])
ax.set_ylim(A_[0],A_[-1])

ax.scatter(neg[1], neg[0], s=0.1, c='white')
if len(points_cvasi) != 0:
    ax.scatter(cva[1], cva[0], s=0.05, c='blue')
ax.scatter(ph[1], ph[0], s=0.1, c='orange')
if len(points_nan) != 0:
    ax.scatter(nan[1], nan[0], s=0.1, c='gray')
plt.show()