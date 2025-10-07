import matplotlib.pyplot as plt
import numpy as np

fig, ax1 = plt.subplots()

x = np.linspace(0, 100)
y1 = np.sin(x)
y2 = np.sqrt(x)

ax1.plot(x, y1)
ax1.set_xlabel("x (основная)")

ax2 = ax1.twiny()
ax2.set_xlabel("x (вторая)")
ax2.plot(x, y2)  # накладывается по x

plt.show()