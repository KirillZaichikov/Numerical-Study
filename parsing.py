import json
import math

def load_end_data(filename):
    """Считывает JSON и парсит данные под ключом 'end_data'."""
    with open(filename, 'r') as f:
        data = json.load(f)

    raw_lines = data["end_data"].strip().split("\n")

    parsed = []
    for line in raw_lines:
        nums = list(map(float, line.split()))
        # nums[0] и nums[1] — параметры
        parsed.append(nums)

    return parsed


def find_closest(parsed_data, p1_target, p2_target):
    """Находит запись с параметрами (p1, p2), ближайшими к заданным."""
    best = None
    best_dist = float("inf")

    for row in parsed_data:
        p1, p2 = row[0], row[1]
        dist = math.hypot(p1 - p1_target, p2 - p2_target)

        if dist < best_dist:
            best_dist = dist
            best = row

    return best


# ====== пример использования ======

filename = "Lyapunov_CelticStoneWN_6D_flow_14-10-2025_23-16-51.json"
p1_target = 0.465       # сюда вставьте значение вашего параметра
p2_target = 735.6      # сюда вставьте второе значение

parsed = load_end_data(filename)
closest = find_closest(parsed, p1_target, p2_target)

print("Ближайшая запись:", closest)
print("Её параметры: p1 =", closest[0], ", p2 =", closest[1])
