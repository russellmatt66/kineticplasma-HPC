# Visualize the particle shapes for various weighting-orders
# Higher-order shape functions are convolutions of lower-order shapes, but np.convolve doesn't seem to be cooperating
import numpy as np
import math
import matplotlib.pyplot as plt


x_grid = np.linspace(-np.pi,np.pi,num=250)

a_0 = (x_grid[x_grid.size - 1] - x_grid[0]) / x_grid.size

shapeFig, shapeAx = plt.subplots(nrows=1, ncols=2)

# 0th-order
h_sq = 1.0 / a_0
S_0 = np.array([h_sq if math.fabs(xj) < a_0 / 2.0 else 0.0 for xj in x_grid])

shapeAx[0].plot(x_grid,S_0)

# 1st-order
h_tri = 2.0 / a_0
S_1 = np.convolve(S_0,S_0,mode='same')
tri_wave = np.empty(x_grid.size)

for j in range(tri_wave.size):
    if (x_grid[j] <= 0 and x_grid[j] > -a_0):
        value = (h_tri / a_0) * (x_grid[j] + a_0)
    elif (x_grid[j] > 0 and x_grid[j] <= a_0):
        value = (h_tri / a_0) * (-x_grid[j] + a_0)
    else:
        value = 0
    tri_wave[j] = value
# print(tri_wave)

# shapeAx[1].plot(x_grid,tri_wave)
shapeAx[1].plot(x_grid,S_1)

# Higher-order 
# S_n(x) = convolve(S_{n-1}, S_{n-1})
S_2 = np.convolve(S_1,S_1,mode='same')

plt.show()