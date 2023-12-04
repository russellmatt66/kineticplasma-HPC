"""
Visualize the particle shapes for various weighting-orders.

This is just a visualization, in a real simulation the diameter of the particle, a_0 = \Delta x, is set equal to the grid-spacing.
Here, that is very difficult to implement, as the 0th-order shape: a square-wave, has issues being constructed for such a fine spacing.
"""
import numpy as np
import math
import matplotlib.pyplot as plt


x_grid = np.linspace(-np.pi,np.pi,num=250)
# print(x_grid)

# a_0 = 1.01 * (x_grid[x_grid.size - 1] - x_grid[0]) / x_grid.size
a_0 = 0.1
print(a_0)

shapeFig, shapeAx = plt.subplots(nrows=2, ncols=2)

# 0th-order
h_sq = 1.0 / a_0
print(h_sq)

S_0 = np.array([h_sq if math.fabs(xj) < a_0 / 2.0 else 0.0 for xj in x_grid])
# print(S_0)

# 1st-order
S_1 = np.convolve(S_0,S_0,mode='same')
# print(S_1)

# tri_wave = np.empty(x_grid.size)
# h_tri = 2.0 / a_0
# for j in range(tri_wave.size):
#     if (x_grid[j] <= 0 and x_grid[j] > -a_0):
#         value = (h_tri / a_0) * (x_grid[j] + a_0)
#     elif (x_grid[j] > 0 and x_grid[j] <= a_0):
#         value = (h_tri / a_0) * (-x_grid[j] + a_0)
#     else:
#         value = 0
#     tri_wave[j] = value
# # print(tri_wave)

shapeAx[0,0].plot(x_grid,S_0)
shapeAx[0,0].set_title('$S_{0}$')
shapeAx[0,0].set_xlim([x_grid[0], x_grid[x_grid.size - 1]])
shapeAx[0,0].set_xlabel('x')
shapeAx[0,0].set_ylabel('A')

# shapeAx[1].plot(x_grid,tri_wave)
shapeAx[0,1].plot(x_grid,S_1)
shapeAx[0,1].set_title('$S_{1}$')
shapeAx[0,1].set_xlim([x_grid[0], x_grid[x_grid.size - 1]])
shapeAx[0,1].set_xlabel('x')
shapeAx[0,1].set_ylabel('A')

S_2 = np.convolve(S_1,S_1,mode='same')
shapeAx[1,0].plot(x_grid,S_2 / max(S_2))
shapeAx[1,0].set_title('$S_{2}$')
shapeAx[1,0].set_xlim([x_grid[0], x_grid[x_grid.size - 1]])
shapeAx[1,0].set_xlabel('x')
shapeAx[1,0].set_ylabel('A')

S_3 = np.convolve(S_2,S_2,mode='same')
shapeAx[1,1].plot(x_grid,S_3 / max(S_3))
shapeAx[1,1].set_title('$S_{3}$')
shapeAx[1,1].set_xlim([x_grid[0], x_grid[x_grid.size - 1]])
shapeAx[1,1].set_xlabel('x')
shapeAx[1,1].set_ylabel('A')

shapeFig.suptitle('Particle shape functions, $S_{m}(x)$, for $m \in [0,3] \cup \mathbb{Z}$')

# Higher-order 
# S_n(x) = convolve(S_{n-1}, S_{n-1})
S_4 = np.convolve(S_3,S_3,mode='same')
S_5 = np.convolve(S_4,S_4,mode='same')
S_6 = np.convolve(S_5,S_5,mode='same')
allFig, allAx = plt.subplots()
allAx.plot(x_grid,S_0 / max(S_0),label='$S_0$')
allAx.plot(x_grid,S_1 / max(S_1),label='$S_1$')
allAx.plot(x_grid,S_2 / max(S_2),label='$S_2$')
allAx.plot(x_grid,S_3 / max(S_3),label='$S_3$')
allAx.plot(x_grid,S_4 / max(S_4),label='$S_4$')
allAx.plot(x_grid,S_5 / max(S_5),label='$S_5$')
allAx.plot(x_grid,S_6 / max(S_6),label='$S_6$')

allAx.set_title('Particle shape-functions, $S_{m}(x)$, for $m \in [0,6] \cup \mathbb{Z}$')
allAx.set_xlabel('x')
allAx.set_ylabel('A')
allAx.set_xlim([x_grid[0], x_grid[x_grid.size - 1]])

plt.legend()
plt.show()