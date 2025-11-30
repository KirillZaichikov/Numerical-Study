import matplotlib.pyplot as plt
import numpy as np
import math as m

A0 = -1.73818
mu0 = -0.58034
alpha = -0.698132
mu_ = np.linspace(1.58, 1.601, 1000)#0.104
A_ = np.linspace(-0.2, 0.09, 1000)
print(mu_+A_)
# new_mu = m.cos(alpha)*(mu_-mu0) - m.sin(alpha)*(A_-A0) + mu0
# new_A = m.sin(alpha)*(mu_-mu0) + m.cos(alpha)*(A_-A0) + A0
# new_mu = 0.76604444311*(mu_-mu0) + 0.64278760968*(A_-A0) + mu0
# new_A = -0.64278760968*(mu_-mu0) + 0.76604444311*(A_-A0) + A0
# print(new_mu)
# mu_ = new_mu
# A_ = new_A


mu_ = np.linspace(-2, 0, 1500)#0.104
A_ = np.linspace(-4, 0, 1500)
C_ = 1.3
nu_ = 0.8#0.964
skip_time = 0
time = 100

times = 1
skip_time = 0
x_ = 0.0
lx, ly = [], []

def map(x, A, C, nu, mu):
    # return abs(-mu+A*x**nu+C*x**(2*nu))
    # return  abs( - ( 0.76604444311 * ( mu + 0.58034 ) + 0.64278760968 * ( A + 1.73818 ) - 0.58034 ) + ( ( - 0.64278760968 ) * ( mu + 0.58034 ) + 0.76604444311 * ( A + 1.73818 ) - 1.73818 ) * x ** nu + C * x ** ( 2.0 * nu ) )
    return abs( - ( 0.5 * ( mu - 1.598 ) - 0.86602540378 * ( A ) + 1.598 ) + ( ( 0.86602540378 ) * ( mu - 1.598 ) + 0.5 * ( A ) ) * x ** nu + C * x ** ( 2.0 * nu ) )

def map_der(x, A, C, nu, mu):
    if map(x, A, C, nu, mu)>0:
        return ( 2.0 * C * nu * pow ( x , 2.0 * nu ) / x + nu * pow ( x , nu ) * ( 0.5 * A + 0.86602540378 * mu - 1.38390859524044 ) / x ) 
        # return ( 2.0 * C * nu * pow ( x , 2.0 * nu ) / x - nu * pow ( x , nu ) * ( - 0.76604444311 * A + 0.64278760968 * mu + 0.779692231276752 ) / x )
    else:
        # return -( 2.0 * C * nu * pow ( x , 2.0 * nu ) / x - nu * pow ( x , nu ) * ( - 0.76604444311 * A + 0.64278760968 * mu + 0.779692231276752 ) / x )
        return -( 2.0 * C * nu * pow ( x , 2.0 * nu ) / x + nu * pow ( x , nu ) * ( 0.5 * A + 0.86602540378 * mu - 1.38390859524044 ) / x ) 
    # if -mu+A*x**nu+C*x**(2*nu) > 0:
        # return A*nu * x **(nu-1) + 2 * nu * C * x ** (2*nu-1)
    # else:
        # return -(A*nu * x **(nu-1) + 2 * nu * C * x ** (2*nu-1))
    
points_pos = []
points_neg = []
points_nan = []
points_cvasi = []
for k1, A in enumerate(A_):
    print('one line is ready')
    for k2, mu in enumerate(mu_):
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
        if min_deriv > 0.008:
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

fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlim(mu_[0],mu_[-1])
ax.set_ylim(A_[0],A_[-1])
# ax.set_xlim(mu_[0],mu_[-1])
# ax.set_ylim(A_[0],A_[-1])

# ax.scatter(pos[1], pos[0], s=0.1, c='white')
# ax.scatter(cva[1], cva[0], s=0.05, c='blue')
ax.scatter(neg[1], neg[0], s=0.1, c='black')
if len(points_nan) != 0:
    ax.scatter(nan[1], nan[0], s=0.1, c='gray')
plt.show()