import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

A = np.float64(0)
nu = np.float64(1)
C = np.float64(1)
np.set_printoptions(precision=17)
def f(x, mu):
    return np.float64(mu) * np.float64(x) * (np.float64(1) - np.float64(x))

def df(x, mu):
    return np.float64(mu) - np.float64(2) * np.float64(mu) * np.float64(x)

def iterate_map(f, x0, mu, n):
    x = np.float64(x0)
    for _ in range(n):
        x = f(x, mu)
    return np.float64(x)

def iterate_list(f, x0, mu, n_transient, n_keep):
    x = np.float64(x0)
    for _ in range(n_transient):
        x = f(x, mu)

    arr = []
    for _ in range(n_keep):
        x = f(x, mu)
        arr.append(x)

    return np.array(arr, dtype=np.float64), np.float64(x)

def detect_period(f, mu, x0, n_transient, period_i_search, tol):
    orbit, x = iterate_list(f, x0, mu, n_transient, period_i_search)
    print(orbit)

    n = len(orbit)
    for p in range(1, n + 1):
        ok = True

        for i in range(n - p):
            if np.abs(orbit[i] - orbit[i + p]) > np.float64(tol):
                print(np.abs(orbit[i] - orbit[i + p]))
                ok = False
                break

        if ok:
            return p, x

    return None, None

def iterate_derivative(f, x0, n):
    x = np.float64(x0)
    deriv = np.float64(1.0)

    for _ in range(n):
        deriv *= df(x)
        x = f(x)

    return np.float64(deriv)

def find_bifurcation(f, mu_left, mu_right, tol_bifurcation, x0, n_transient, period_i_search, tol_detect_period):
    for _ in range(100):
        mu_mid = np.float64(0.5) * (np.float64(1.4) * mu_left + np.float64(0.6) * mu_right)

        p, x0 = detect_period(f, mu_mid, x0, n_transient, period_i_search, tol_detect_period)

        if p < period_i_search:
            mu_left = mu_mid
        elif p == period_i_search:
            print("test", p, mu_mid)
            mu_right = mu_mid
        else:
            print("ERROR")
            break

        if np.abs(mu_right - mu_left) < np.float64(tol_bifurcation):
            print(mu_left, mu_right)
            print(_)
            break

    return np.float64(0.5) * (mu_left + mu_right)

def feigenbaum_cascade(f, mu0, x00, skip_iter_num, tol_det_per, tol_bif, search_width=0.1, levels=8):
    mus = [np.float64(mu0)]
    current_mu = np.float64(mu0)
    target_period = 1

    for n in range(levels):
        print(f"\nSearching bifurcation for period {target_period} -> {2*target_period}")

        mu_left = current_mu
        mu_right = current_mu + np.float64(search_width)

        i = 0
        while True:
            i += 1
            p_right, x00 = detect_period(
                f, mu_right, x00, skip_iter_num, 2 * target_period, tol_det_per
            )

            print(p_right)

            if p_right is not None and p_right == 2 * target_period:
                break
            elif p_right is not None and p_right > 2 * target_period:
                mu_right -= np.float64(search_width)
                search_width /= np.float64(2)
                mu_right += np.float64(search_width)
            elif p_right is not None and p_right < 2 * target_period:
                mu_right += np.float64(search_width)

        mu_new = find_bifurcation(
            f,
            mu_left,
            mu_right,
            tol_bif,
            x00,
            skip_iter_num,
            2 * target_period,
            tol_det_per
        )

        mus.append(mu_new)

        current_mu = mu_new
        target_period *= 2

        print(f"mu_{len(mus)-1} = {mu_new:.15f}")
        break

    return np.array(mus, dtype=np.float64)

def estimate_delta(mus):
    mus = np.array(mus, dtype=np.float64)

    deltas = []
    for n in range(2, len(mus)-1):
        d = (mus[n] - mus[n-1]) / (mus[n+1] - mus[n])
        deltas.append(d)

    return np.array(deltas, dtype=np.float64)


if __name__ == "__main__":
    mu0 = np.float64(1.5)
    x00 = np.float64(0.1)
    skip_time = 1000

    mus = feigenbaum_cascade(
        f,
        mu0=mu0,
        x00=x00,
        skip_iter_num=skip_time,
        tol_det_per=np.float64(1e-8),
        tol_bif=np.float64(1e-12),
        search_width=np.float64(0.05),
        levels=7
    )

    print("\nBifurcation points:")
    for i, mu in enumerate(mus):
        print(f"mu_{i} = {mu:.15f}")

    deltas = estimate_delta(mus)

    print("\nFeigenbaum delta estimates:")
    for i, d in enumerate(deltas):
        print(f"delta_{i} = {d:.12f}")

    print("\nExpected:")
    print("delta = 4.669201609102...")

    plt.plot(range(len(deltas)), deltas, 'o-')
    plt.axhline(4.669201609102, linestyle='--')
    plt.xlabel("n")
    plt.ylabel("delta_n")
    plt.title("Feigenbaum delta convergence")
    plt.grid(True)
    plt.show()