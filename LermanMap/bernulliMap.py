import numpy as np
import math as m
import matplotlib.pyplot as plt
import time
from params import *
from Step import *
from model import *



fig = plt.figure(figsize=(5, 5))

ax2d_2 = fig.add_subplot()
ax2d_2.tick_params(axis='both', labelsize=15)
ticks = [
    -2*np.pi, -3*np.pi/2, -np.pi, -np.pi/2, 
    0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi
]
labels = [
    r'$-2\pi$', r'$-\frac{3\pi}{2}$', r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$',
    r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'
]
ax2d_2.set_xticks(ticks)
ax2d_2.set_xticklabels(labels, fontsize=16)

ax2d_2.set_yticks(ticks)
ax2d_2.set_yticklabels(labels, fontsize=16)

ax2d_2.set_xlim(-np.pi, np.pi)
ax2d_2.set_ylim(-np.pi, np.pi)

def drawing(phi_n, phi_n1, clr):
    # ax2d_2.plot(phi_n, phi_n1, linestyle="", marker="o", markersize=0.05, color=clr, rasterized=True)
    # ax2d_2.plot(phi_n, phi_n1, '.', color=clr, alpha=0.25, markersize=0.03)
    ax2d_2.plot(phi_n, phi_n1, ',', color=clr, alpha=0.25)
    plt.savefig('high_res_plot.png', dpi=400)


def main(start_point, clr):
    n = int((integrate_time / step) / skip_for_phase)
    phi_n = np.zeros(n, dtype=np.float64)
    phi_n1 = np.zeros(n, dtype=np.float64)
    
    print("***********TIME for integrate:", (int(integrate_time / step)))

    start = time.time()
    for i in range(int(skip_time / step)):
        makeStep(start_point, dimension, diffFunc, params, step)
    print("***********END warm time:", time.time() - start)

    start = time.time()
    for i in range(int(integrate_time / step)):
        if i % skip_for_phase == 0:
            phi_n[i // skip_for_phase] = start_point[2]
        makeStep(start_point, dimension, diffFunc, params, step)
        if i % skip_for_phase == 0:
            phi_n1[i // skip_for_phase] = start_point[2]
    end = time.time()
    print("***********TIME:", end-start)

    drawing(phi_n, phi_n1, clr)


if __name__ == "__main__":
    start_point = np.array(initial_point)
    clrs = ["black", "red"]
    print(start_point)
    for (i, j) in zip(start_point, clrs):
        print(i)
        main(i, j)
    print("*********** calc is ending")
    plt.savefig("plot.pdf", dpi=300)
    plt.show()