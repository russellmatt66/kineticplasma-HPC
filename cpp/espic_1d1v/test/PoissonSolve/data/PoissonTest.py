import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

""" 
Parse result datafiles with pandas and plot them along with the analytical solution
"""

# Refactored code test
# Case 1: rho(x) = sin(x) => phi(x) = sin(x)
sin_df = pd.read_csv("sin.csv")
x_grid = sin_df['x_j']

sin_df.plot(x='x_j', y='phi_x', label='Numerical-Refactored')
plt.plot(sin_df['x_j'], sin_df['phi_original'], label='Numerical-Original')
plt.plot(sin_df['x_j'], np.sin(x_grid), label='Analytical: sin(x)')
plt.xlabel('x')
plt.ylabel('phi(x)')
plt.legend()

# Case 2: rho(x) = cos(x) => phi(x) = cos(x)
cos_df = pd.read_csv("cos.csv")

cos_df.plot(x='x_j', y = 'phi_x', label='Numerical-Refactored')
plt.plot(cos_df['x_j'], cos_df['phi_original'], label='Numerical-Original')
plt.plot(cos_df['x_j'], np.cos(x_grid), label='Analytical: cos(x)')
plt.xlabel('x')
plt.ylabel('phi(x)')
plt.legend()

# Case 3: rho(x) = x^2 => phi(x) = -x^4 / 12
para_df = pd.read_csv("parabola.csv")

para_df.plot(x='x_j', y = 'phi_x', label='Numerical-Refactored')
plt.plot(para_df['x_j'], para_df['phi_original'], label='Numerical-Original')
plt.plot(para_df['x_j'], -1.0 / 12.0 * x_grid**4, label='Analytical: $-x^{4} / 12$')
plt.xlabel('x')
plt.ylabel('phi(x)')
plt.legend()

# Case 4: rho(x) = e^x => phi(x) = -e^x
exp_df = pd.read_csv("exponential.csv")

exp_df.plot(x='x_j', y = 'phi_x', label='Numerical-Refactored')
plt.plot(exp_df['x_j'], exp_df['phi_original'], label='Numerical-Original')
plt.plot(exp_df['x_j'], -np.exp(x_grid), label='Analytical: $-e^{x}$')
plt.xlabel('x')
plt.ylabel('phi(x)')
plt.legend()

plt.show()