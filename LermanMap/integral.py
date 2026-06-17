from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

# Вариант 1: В лоб
# p1 = -pow(np.e,  alpha * t) * np.cos(beta * t + teta) * pow(0.2e1, 0.3e1 / 0.4e1) * pow(np.cosh(0.4e1 * alpha * t), -0.1e1 / 0.2e1) / 0.2e1
# p2 = -pow(np.e,  alpha * t) * np.sin(beta * t + teta) * pow(0.2e1, 0.3e1 / 0.4e1) * pow(np.cosh(0.4e1 * alpha * t), -0.1e1 / 0.2e1) / 0.2e1
# q1 =  pow(np.e, -alpha * t) * np.cos(beta * t + teta) * pow(0.2e1, 0.3e1 / 0.4e1) * pow(np.cosh(0.4e1 * alpha * t), -0.1e1 / 0.2e1) / 0.2e1
# q2 =  pow(np.e, -alpha * t) * np.sin(beta * t + teta) * pow(0.2e1, 0.3e1 / 0.4e1) * pow(np.cosh(0.4e1 * alpha * t), -0.1e1 / 0.2e1) / 0.2e1

# Вариант 2: экспоненты с угасающими степениями (возникает divide by zero)
# if t>0:
#     p1 = -np.cos(beta * t + teta) * pow(2, -0.25) * pow(((1 + np.exp(-8*alpha*t)) / 2), -0.5) * (np.exp(-  alpha*t))
#     p2 = -np.sin(beta * t + teta) * pow(2, -0.25) * pow(((1 + np.exp(-8*alpha*t)) / 2), -0.5) * (np.exp(-  alpha*t))
#     q1 =  np.cos(beta * t + teta) * pow(2, -0.25) * pow(((1 + np.exp(-8*alpha*t)) / 2), -0.5) * (np.exp(-3*alpha*t))
#     q2 =  np.sin(beta * t + teta) * pow(2, -0.25) * pow(((1 + np.exp(-8*alpha*t)) / 2), -0.5) * (np.exp(-3*alpha*t))
# else:
#     p1 = -np.cos(beta * t + teta) * pow(2, -0.25) * pow(((1 + 1/np.exp(8*alpha*t)) / 2), -0.5) / (np.exp(  alpha*t))
#     p2 = -np.sin(beta * t + teta) * pow(2, -0.25) * pow(((1 + 1/np.exp(8*alpha*t)) / 2), -0.5) / (np.exp(  alpha*t))
#     q1 =  np.cos(beta * t + teta) * pow(2, -0.25) * pow(((1 + 1/np.exp(8*alpha*t)) / 2), -0.5) / (np.exp(3*alpha*t))
#     q2 =  np.sin(beta * t + teta) * pow(2, -0.25) * pow(((1 + 1/np.exp(8*alpha*t)) / 2), -0.5) / (np.exp(3*alpha*t))

alpha, beta = 1, 1
A = alpha / np.sqrt(2)

def determ(teta):

    def solution(t):
        lnp1 = -  alpha*t + (1/4) * np.log(2) + np.log(abs(np.cos(beta*t + teta))) - (1/2)*np.logaddexp(0, -8*alpha*t)
        lnp2 = -  alpha*t + (1/4) * np.log(2) + np.log(abs(np.sin(beta*t + teta))) - (1/2)*np.logaddexp(0, -8*alpha*t)
        lnq1 = -3*alpha*t + (1/4) * np.log(2) + np.log(abs(np.cos(beta*t + teta))) - (1/2)*np.logaddexp(0, -8*alpha*t)
        lnq2 = -3*alpha*t + (1/4) * np.log(2) + np.log(abs(np.sin(beta*t + teta))) - (1/2)*np.logaddexp(0, -8*alpha*t)
        p1 = - np.sign(np.cos(beta*t + teta)) * np.exp(lnp1)
        p2 = - np.sign(np.sin(beta*t + teta)) * np.exp(lnp2)
        q1 =   np.sign(np.cos(beta*t + teta)) * np.exp(lnq1)
        q2 =   np.sign(np.sin(beta*t + teta)) * np.exp(lnq2)
        return p1, p2, q1, q2

    def func11(t):
        p1, p2, q1, q2 = solution(t)
        # f11 = 2 * (alpha * q1 - beta * q2 + 4 * A * (p1 * p1 + p2 * p2) * p1) * p1 + 2 * (alpha * p1 + beta * p2 + 4 * A * (q1 * q1 + q2 * q2) * q1) * q1
        f11 =  -2 * (alpha * q1 - beta * q2 + 4 * A * (p1 * p1 + p2 * p2) * p1) * q1 + 2 * (alpha * p1 + beta * p2 + 4 * A * (q1 * q1 + q2 * q2) * q1) * p1
        return f11
    
    def func12(t):
        p1, p2, q1, q2 = solution(t)
        f12 = (alpha * q1 - beta * q2 + 4 * A * (p1 * p1 + p2 * p2) * p1) * (-p1 - p2) + (alpha * q2 + beta * q1 + 4 * A * (p1 * p1 + p2 * p2) * p2) * (p1 - p2) + (alpha * p1 + beta * p2 + 4 * A * (q1 * q1 + q2 * q2) * q1) * (-q1 + q2) + (alpha * p2 - beta * p1 + 4 * A * (q1 * q1 + q2 * q2) * q2) * (-q1 - q2)
        return f12
    
    def func21(t):
        p1, p2, q1, q2 = solution(t)
        # f21 = 2 * q2 * p1 - 2 * p2 * q1
        f21 =  -2 * p2 * p1 - 2 * q2 * q1
        return f21

    def func22(t):
        p1, p2, q1, q2 = solution(t)
        f22 = q2 * (-p1 - p2) - q1 * (p1 - p2) - p2 * (-q1 + q2) + p1 * (-q1 - q2)
        return f22
    
    d11, _ = quad(func11, -np.inf, np.inf)
    d12, _ = quad(func12, -np.inf, np.inf)
    d21, _ = quad(func21, -np.inf, np.inf)
    d22, _ = quad(func22, -np.inf, np.inf)
    print(_)
    # print(d11,d12,d21,d22)

    return d11 * d22 - d21 * d12


tetaMas = np.linspace(0,2*np.pi,200)
x = []
y = []
for i, teta0 in enumerate(tetaMas):
    x.append(teta0)
    print(i)
    y.append(determ(teta0))

fig,ax = plt.subplots()
ax.set_xlabel("teta")
ax.plot(x, y)
ax.plot(x, [0]*len(x), color="red")
plt.show()