import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


eps = 1e-3
skip_time = 1000 
DISCR = 1000 # disc param
MAX_PERIOD = 40

# for diagram
fig, ax = plt.subplots()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1, 4)


def map(initial_point, params):
    return params[0] + params[1] * initial_point - initial_point ** 3


def gen_hex(number):
    res = hex(number)[2:]
    if len(res) < 6:
        res = '0' * (6 - len(res)) + res
    elif len(res) > 6:
        res = res[len(res) - 6:]
    return res


color_values = [i/MAX_PERIOD for i in range(MAX_PERIOD+1)]
color_names1 = ['skyblue', 'aqua', 'green', 'yellow', 'red', 'darkblue', 'pink', 
               'grey', 'olive', 'peru', 'moccasin', 'lime', 'teal', 'dimgrey',
               'midnightblue', 'yellowgreen', 'coral', 'maroon', 'orchid', 'cornsilk', 'gold']
color_names2 = [f'#{gen_hex(i*5000)}' for i in range(10, 10 + MAX_PERIOD - 20)]
color_names = color_names1 + color_names2
color_names[-1] = 'black'
cMap = list(zip(color_values, color_names))
print(cMap)
print(cMap)
customColourMap = LinearSegmentedColormap.from_list("pri_c", cMap)
plt.colormaps.register(customColourMap)


def out_bar(percent, first):
    if not first:
        print ("\033[A                             \033[A")
    percent = int(percent)
    text = '[' + percent * '*' + (100 - percent) * '-' + ']' + f'  {percent}%'
    print(text)
    print("\033[33m".format(text), end='')


def main():
    a = np.linspace(-2.2, 2.2, DISCR)
    b = np.linspace(-1.5, 5, DISCR)

    mas_for_points = np.ndarray((DISCR*DISCR, 3))
    old_percent = 0
    first = True
    print('calc start:')
    for i in range(DISCR):
        # initial_point = eps
        for j in range(DISCR):
            params = [a[i], b[j]]
            initial_point = 0

            for counter in range(skip_time):
                initial_point = map(initial_point, params)

            first_point = initial_point
            k = 0
            while k < MAX_PERIOD:
                k += 1
                initial_point = map(initial_point, params)
                if (initial_point <= first_point + eps) and (initial_point >= first_point - eps):
                    break
            mas_for_points[i * DISCR + j] = [params[0], params[1], k / MAX_PERIOD]

            current_percent = (100 * ((i*DISCR) + j) / (DISCR * DISCR)) // 1
            if current_percent > old_percent:
                out_bar((100 * ((i*DISCR) + j) / (DISCR * DISCR)) // 1, first)
                first = False
            old_percent = current_percent

    mas_for_points = mas_for_points[~np.isnan(mas_for_points).any(axis=1)]
    mas_for_points = mas_for_points.T

    ax.scatter(mas_for_points[0], mas_for_points[1], s=1, c=mas_for_points[2], cmap='pri_c') # prism
    plt.show()


if __name__ == "__main__":
    main()