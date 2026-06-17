import numpy as np

def f(x, r):
    return r * x * (1 - x)

def dfdx(x, r, h=1e-8):
    return (f(x+h, r) - f(x-h, r)) / (2*h)

def iterate(x, r, k):
    for _ in range(k):
        x = f(x, r)
    return x

def diterate(x, r, k):
    p = 1.0
    for _ in range(k):
        p *= dfdx(x, r)
        x = f(x, r)
    return p

def newton(r, k, x, tol=1e-14):
    for _ in range(100):
        fx = iterate(x, r, k) - x
        dfx = diterate(x, r, k) - 1
        if abs(dfx) < 1e-14:
            return None
        dx = fx / dfx
        x -= dx
        if abs(dx) < tol:
            return x
    return None

def bifurcation(r, k, x):
    x = newton(r, k, x)
    return None if x is None else diterate(x, r, k) + 1

def find_r(k, r, x, dr=1e-5):
    g0 = bifurcation(r, k, x)
    while True:
        r2 = r + dr
        g1 = bifurcation(r2, k, x)
        if g1 is None:
            r = r2
            continue
        if g0 * g1 < 0:
            break
        r, g0 = r2, g1

    a, b = r, r2
    for _ in range(60):
        m = 0.5 * (a + b)
        gm = bifurcation(m, k, x)
        if gm is None:
            continue
        if abs(gm) < 1e-14:
            return m
        if bifurcation(a, k, x) * gm < 0:
            b = m
        else:
            a = m
    return 0.5 * (a + b)

r0 = 2.5
x0 = 0.5
k = 1
rs = []

for _ in range(8):
    r1 = find_r(k, r0, x0)
    rs.append(r1)
    print(k, r1)
    x0 = newton(r1, 2*k, x0)
    r0 = r1 + 1e-6
    k *= 2

for i in range(len(rs)-2):
    print((rs[i+1]-rs[i])/(rs[i+2]-rs[i+1]))