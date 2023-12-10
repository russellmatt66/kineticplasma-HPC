import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

""" 
Parse result datafiles with pandas and plot them along with the analytical solution
"""

""" Violating DRY, should be refactored with a wrapper for the plotting """

# Refactored code test
# Case 1: rho(x) = sin(x) => phi(x) = sin(x)
sin_df = pd.read_csv("sin.csv")
x_grid = sin_df['x_j']

sinFig, sinAx = plt.subplots(nrows=1,ncols=2)
sinFig.suptitle('Case 1: rho(x) = sin(x)')
# phi(x)
sinAx[0].plot(x_grid, sin_df['phi_x'], label='Numerical-Refactored', color='blue')
sinAx[0].plot(x_grid, sin_df['phi_original'], '--', label='Numerical-Original', color='orange')
sinAx[0].plot(x_grid, np.sin(x_grid), label='Analytical: $sin(x)$', color='green')
sinAx[0].set_xlabel('x')
sinAx[0].set_ylabel('phi(x)')
sinAx[0].legend()

# E(x)
sinAx[1].plot(x_grid, sin_df['E_j'], label='Numerical-Refactored', color='blue')
sinAx[1].plot(x_grid, -np.cos(x_grid), label='Analytical: $-cos(x)$', color='green')
sinAx[1].set_xlabel('x')
sinAx[1].set_ylabel('E(x)')
sinAx[1].legend()

# Case 2: rho(x) = cos(x) => phi(x) = cos(x)
cos_df = pd.read_csv("cos.csv")

cosFig, cosAx = plt.subplots(nrows=1,ncols=2)
cosFig.suptitle('Case 2: rho(x) = cos(x)')
# phi(x)
cosAx[0].plot(x_grid, cos_df['phi_x'], label='Numerical-Refactored', color='blue')
cosAx[0].plot(x_grid, cos_df['phi_original'], '--', label='Numerical-Original', color='orange')
cosAx[0].plot(x_grid, np.cos(x_grid), label='Analytical: $cos(x)$', color='green')
cosAx[0].set_xlabel('x')
cosAx[0].set_ylabel('phi(x)')
cosAx[0].legend()

# E(x)
cosAx[1].plot(x_grid, cos_df['E_j'], label='Numerical-Refactored', color='blue')
cosAx[1].plot(x_grid, np.sin(x_grid), label='Analytical: $sin(x)$', color='green')
cosAx[1].set_xlabel('x')
cosAx[1].set_ylabel('E(x)')
cosAx[1].legend()


# Case 3: rho(x) = x^2 => phi(x) = -x^4 / 12
para_df = pd.read_csv("parabola.csv")

paraFig, paraAx = plt.subplots(nrows=1,ncols=2)
paraFig.suptitle('Case 3: rho(x) = $x^{2}$')
# phi(x)
paraAx[0].plot(x_grid, para_df['phi_x'], label='Numerical-Refactored', color='blue')
paraAx[0].plot(x_grid, para_df['phi_original'], '--', label='Numerical-Original', color='orange')
paraAx[0].plot(x_grid, -1.0 / 12.0 * x_grid**4, label='Analytical: $-x^{4} / 12$', color='green')
paraAx[0].set_xlabel('x')
paraAx[0].set_ylabel('phi(x)')
paraAx[0].legend()

# E(x)
paraAx[1].plot(x_grid, para_df['E_j'], label='Numerical-Refactored', color='blue')
paraAx[1].plot(x_grid, x_grid**3 / 3.0, label='Analytical: $x^{3} / 3$', color='green')
paraAx[1].set_xlabel('x')
paraAx[1].set_ylabel('E(x)')
paraAx[1].legend()

# Case 4: rho(x) = e^x => phi(x) = -e^x
exp_df = pd.read_csv("exponential.csv")

expFig, expAx = plt.subplots(nrows=1,ncols=2)
expFig.suptitle('Case 4: rho(x) = exp(x)')
# phi(x)
expAx[0].plot(x_grid, exp_df['phi_x'], label='Numerical-Refactored')
expAx[0].plot(x_grid, exp_df['phi_original'], '--', label='Numerical-Original')
expAx[0].plot(x_grid, -np.exp(x_grid), label='Analytical: $-e^{x}$')
expAx[0].set_xlabel('x')
expAx[0].set_ylabel('phi(x)')
expAx[0].legend()

# E(x)
expAx[1].plot(x_grid, exp_df['E_j'], label='Numerical-Refactored')
expAx[1].plot(x_grid, np.exp(x_grid), label='Analytical: $e^{x}$')
expAx[1].set_xlabel('x')
expAx[1].set_ylabel('E(x)')
expAx[1].legend()

plt.show()