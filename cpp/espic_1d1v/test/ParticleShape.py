# Visualize the particle shapes for various weighting-orders
import numpy as np
import math
import matplotlib.pyplot as plt

a_0 = 2.0 # particle diameter
x_grid = np.linspace(-2.0*np.pi,2.0*np.pi,num=250)

shapeFig, shapeAx = plt.subplots(nrows=1, ncols=2)

# 0th-order
h_sq = 1.0 / a_0
square_wave = np.array([h_sq if math.fabs(xj) < a_0 / 2.0 else 0.0 for xj in x_grid])
# print(square_wave)

shapeAx[0].plot(x_grid,square_wave)

# 1st-order
h_tri = 0.5 * (1.0 + a_0 / 2.0)
tri_wave = np.empty(x_grid.size)

for j in range(tri_wave.size):
    if (x_grid[j] < 0 and x_grid[j] >= -a_0 / 2.0):
        value = (2.0 / a_0) * (x_grid[j] + h_tri)
    elif (x_grid[j] > 0 and x_grid[j] <= a_0 / 2.0):
        value = (2.0 / a_0) * (-x_grid[j] + h_tri)
    else:
        value = 0
    tri_wave[j] = value
# print(tri_wave)

shapeAx[1].plot(x_grid,tri_wave)

# Higher-order

plt.show()