import json
import numpy as np
import matplotlib.pyplot as plt

TARGET_P2 = 0.05


def safe_float(x):
    """
    Преобразует строку в float.
    Любые невалидные значения (-nan(ind), nan(ind), 1.#IND и т.п.)
    превращаются в np.nan.
    """
    try:
        return float(x)
    except ValueError:
        return np.nan


def kaplan_yorke_dimension(exponents):
    """
    Размерность Каплана-Йорке.
    """
    lam = np.sort(np.asarray(exponents))[::-1]

    cumulative = 0.0

    for j in range(len(lam)):
        cumulative += lam[j]

        if cumulative < 0:
            if j == 0:
                return 0.0

            s_prev = cumulative - lam[j]
            return j + s_prev / abs(lam[j])

    return float(len(lam))


# ---------- чтение JSON ----------

with open("Lyapunov_lerFRMS2_3D_map_17-06-2026_16-19-49.json", "r", encoding="utf-8") as f:
    data = json.load(f)

rows = []

for line in data["_Data"].splitlines():
    line = line.strip()

    if not line:
        continue

    vals = [safe_float(x) for x in line.split()]

    if len(vals) < 5:
        continue

    rows.append(vals[:5])

rows = np.array(rows)

if len(rows) == 0:
    raise RuntimeError("В файле не найдено ни одной корректной строки данных.")

# ---------- удаляем строки с NaN в показателях ----------

valid_mask = ~np.isnan(rows[:, 2:5]).any(axis=1)
rows = rows[valid_mask]

if len(rows) == 0:
    raise RuntimeError(
        "После удаления строк с NaN не осталось ни одной точки."
    )

# ---------- поиск ближайшего p2 ----------

all_p2 = rows[:, 1]
unique_p2 = np.unique(all_p2)

selected_p2 = unique_p2[np.argmin(np.abs(unique_p2 - TARGET_P2))]

print(f"Искомое p2 = {TARGET_P2}")
print(f"Ближайшее значение в сетке = {selected_p2}")

# ---------- выбор горизонтали ----------

mask = rows[:, 1] == selected_p2
section = rows[mask]

if len(section) == 0:
    raise RuntimeError(
        f"Не удалось найти точки для p2 = {selected_p2}"
    )

# сортировка по p1
section = section[np.argsort(section[:, 0])]

# ---------- размерность Каплана-Йорке ----------

p1_values = []
dimensions = []

for row in section:
    lyap = row[2:5]

    D = kaplan_yorke_dimension(lyap)

    if row[0] > 0.073:
        p1_values.append(row[0])
        dimensions.append(D)

p1_values = np.array(p1_values)
dimensions = np.array(dimensions)

# ---------- график ----------

plt.figure(figsize=(10, 6))
plt.plot(p1_values, dimensions, linewidth=3)
plt.axhline(2, c='red', ls='--')
plt.axhline(1, c='red', ls='--')

plt.axvline(0.15, c='black', ls='--')
plt.axvline(0.3, c='black', ls='--')
plt.axvline(0.48, c='black', ls='--')
plt.axvline(0.7, c='black', ls='--')

plt.text(0.16, 1.1, "a", fontsize=18)
plt.text(0.31, 1.1, "b", fontsize=18)
plt.text(0.49, 1.1, "c", fontsize=18)
plt.text(0.71, 1.1, "d", fontsize=18)

ax = plt.gca()
ax.set_xlabel(r"$\alpha$", fontsize=16)
ax.set_ylabel(r"$D_{KY}$", fontsize=16)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.grid(True)

plt.tight_layout()
plt.show()