import matplotlib.pyplot as plt
import numpy as np

LYAP_TRES = 0.001
DERI_TRES = 0.005

mu_ = np.linspace(-0.3, 3, 500)
A_ = np.linspace(-1.5, 1.1, 500)
C_ = 1.3
nu_ = 1
time = 1000
skip_time = 200


A_grid, mu_grid = np.meshgrid(A_, mu_, indexing='ij')

def map_vec(x, A, C, nu, mu):
    return np.abs(-mu + A * x**nu + C * x**(2*nu))

def map_der_vec(x, A, C, nu, mu):
    inner = -mu + A * x**nu + C * x**(2*nu)
    deriv = A * nu * x**(nu - 1) + 2 * nu * C * x**(2*nu - 1)
    return np.where(inner > 0, deriv, -deriv)


x = np.zeros_like(A_grid)

for _ in range(skip_time):
    x = map_vec(x, A_grid, C_, nu_, mu_grid)

lyap = np.zeros_like(A_grid)
min_deriv = np.full_like(A_grid, np.inf)
has_inf = np.zeros_like(A_grid, dtype=bool)

for _ in range(time):
    x = map_vec(x, A_grid, C_, nu_, mu_grid)
    d = map_der_vec(x, A_grid, C_, nu_, mu_grid)
    abs_d = np.abs(d)

    bad = abs_d <= 0
    has_inf |= bad

    safe = ~bad
    lyap = np.where(safe, lyap + np.log(np.where(safe, abs_d, 1.0)), lyap)
    min_deriv = np.minimum(min_deriv, np.where(safe, abs_d, np.inf))

is_nan   = has_inf
is_ph    = ~is_nan & (min_deriv >  DERI_TRES) & (lyap >  LYAP_TRES)
is_cvasi = ~is_nan & (min_deriv <= DERI_TRES) & (lyap >  LYAP_TRES)
is_neg   = ~is_nan & (lyap <= LYAP_TRES)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(mu_[0], mu_[-1])
ax.set_ylim(A_[0], A_[-1])

ax.scatter(mu_grid[is_neg],   A_grid[is_neg],   s=0.1,  c='white')
ax.scatter(mu_grid[is_cvasi], A_grid[is_cvasi], s=0.05, c='blue')
ax.scatter(mu_grid[is_ph],    A_grid[is_ph],    s=0.1,  c='orange')
ax.scatter(mu_grid[is_nan],   A_grid[is_nan],   s=0.1,  c='gray')

ax.set_facecolor('gray')
plt.tight_layout()
plt.show()