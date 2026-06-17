import numpy as np
from scipy.integrate import quad

# Производные координат
def dx_dt(t):
    # return 2
    return 1

def dy_dt(t):
    # return -(1/2) * np.sin(t)
    return (1) * np.sin(t)

def dz_dt(t):
    # return (1/2) * np.cos(t)
    return (1) * np.cos(t)

# Подынтегральная функция
def ds(t):
    return np.sqrt(
        dx_dt(t)**2 +
        dy_dt(t)**2 +
        dz_dt(t)**2
    )

# Пределы интегрирования
t0 = 0
t1 = 2 * np.pi

# Численное интегрирование
L, error = quad(ds, t0, t1)

print("Длина кривой =", L)
print("Оценка ошибки =", error)
print("Длина кривой / 2*np.pi =", L / (2*np.pi))


'''
При дествии отображением смейла вильямса из статьи Кузнецова, 
то длина кривой (0,0,phi) увеличивается примерно в 2 раза


'''