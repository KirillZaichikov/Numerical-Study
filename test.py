
import math as m

A0 = -1.73818
mu0 = -0.58034
alpha = -m.pi/6


def f(A, mu):
    new_mu = m.cos(alpha)*(mu-mu0) - m.sin(alpha)*(A-A0) + mu0
    new_A = m.sin(alpha)*(mu-mu0) + m.cos(alpha)*(A-A0) + A0
    return new_mu, new_A


A = -2.5
mu = -1
print(f(A,mu))