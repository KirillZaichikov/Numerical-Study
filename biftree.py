import numpy as np
import math as m
import matplotlib.pyplot as plt


eps = 1e-3
skip_time = 1000
DISCR = 1000 # disc param
a = 0.1
b = np.linspace(-0.5, 2.8, DISCR)

fig, ax = plt.subplots()


def map(initial_point, params):
    return params[0] + params[1] * initial_point - initial_point ** 3


def main():

    mas_for_points = np.ndarray((DISCR*50, 2))
    print('calc start:')
    initial_point = 0.25  # надо указывать
    for i in range(DISCR):
        params = [a, b[DISCR - i-1]]
        print(params)

        for counter in range(skip_time):
            initial_point = map(initial_point, params)

        for counter in range(50):
            initial_point = map(initial_point, params)
            mas_for_points[counter+i*50] = [params[1], initial_point]

    mas_for_points = mas_for_points[~np.isnan(mas_for_points).any(axis=1)]
    mas_for_points = mas_for_points[~np.isinf(mas_for_points).any(axis=1)]
    print(mas_for_points)
    mas_for_points = mas_for_points.T
    ax.scatter(mas_for_points[0], mas_for_points[1], s=1, c="black")

    mas_for_points = np.ndarray((DISCR*50, 2))
    initial_point = 0.25  # надо указывать
    for i in range(DISCR):
        params = [a, b[i]]
        print(params)

        for counter in range(skip_time):
            initial_point = map(initial_point, params)

        for counter in range(50):
            initial_point = map(initial_point, params)
            mas_for_points[counter+i*50] = [params[1], initial_point]

    mas_for_points = mas_for_points[~np.isnan(mas_for_points).any(axis=1)]
    mas_for_points = mas_for_points[~np.isinf(mas_for_points).any(axis=1)]
    print(mas_for_points)
    mas_for_points = mas_for_points.T
    ax.scatter(mas_for_points[0], mas_for_points[1], s=1, c="skyblue")
    plt.show()


if __name__ == "__main__":
    main()