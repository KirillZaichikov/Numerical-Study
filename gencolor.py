import matplotlib.pyplot as plt
import numpy as np
import distinctipy


colors = distinctipy.get_colors(800, pastel_factor=0.2) # , pastel_factor=0.4
print(colors)
colormap = np.array(colors)
distinctipy.color_swatch(colormap)
plt.show()