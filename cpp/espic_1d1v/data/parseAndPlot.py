import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

from matplotlib.animation import FuncAnimation
"""
Read datafiles in Grid/ and Particles/
"""
grid_path = 'Grid'
particle_path = 'Particles'

# Creates a list of all relevant files based on the argument
grid_files = glob.glob(grid_path + "/*.csv")
particle_files = glob.glob(particle_path + "/*.csv")

# Little slower than using generator expression, i.e., () instead of []
grid_df_list = [pd.read_csv(datafile) for datafile in grid_files]
particle_df_list = [pd.read_csv(datafile) for datafile in particle_files]

# print(particle_df_list.type())

grid_df = pd.concat(grid_df_list, ignore_index=True)
particle_df = pd.concat(particle_df_list, ignore_index=True)

# Looks good here
# print(grid_df)
# print(particle_df)

"""
Parse DataFrames
"""
# Integrate electric field and particle kinetic energy
def trapezoidalIntegration(f: pd.Series, dx: float) -> float:
    sum = 0.0
    for j in np.arange(0,f.size - 1): # This is getting short-circuited
        sum += 0.5 * (f.values[j] + f.values[j+1])
    sum *= dx
    print(sum)
    return sum

# Appending to lists is just the simplest solution I can think of
ElectricEnergy = []
ParticleKineticEnergy = []
TotalE = []

first_grid_df = grid_df_list[0]
dx = (first_grid_df['x_j'].max() - first_grid_df['x_j'].min()) / (first_grid_df['j'].max()) 
print(dx)

x_grid = list(first_grid_df['x_j'])

time_vec = np.arange(len(grid_df_list))

# Probably slow
for it in time_vec:
    current_grid_df = grid_df_list[it]
    current_particle_df = particle_df_list[it]
    ElectricEnergy.append(trapezoidalIntegration(0.5 * current_grid_df['E_j'] * current_grid_df['E_j'],dx)) # <- Core problem here! 'numpy.float64' 
    KE = 0.5 * current_particle_df['v_i'].dot(current_particle_df['v_i']) # 
    ParticleKineticEnergy.append(KE)
    TotalE.append(KE + ElectricEnergy[it])


"""
Plot output
"""
# Create plot of energy history
# - Looks weird for N = 1024, Nx = 2048
# - Takes strange dips
energyFig, energyAx = plt.subplots()

energyAx.plot(time_vec, TotalE, label='TotalE')
energyAx.plot(time_vec, ElectricEnergy, label='E')
energyAx.plot(time_vec, ParticleKineticEnergy, label='KE')
energyAx.legend()

# Create movie of grid field and potential

# Create movie of phase-space

# Create movie of distribution function
# Need to make histogram into distribution function
# Make f(v) out of histogram by normalizing counts
"""
fig, ax = plt.subplots()

v_min = particle_df_list[0]['v_i'].min()
v_max = particle_df_list[0]['v_i'].max()

def update(frame):
    ax.clear()
    particle_df_frame = particle_df_list[frame]
    N = particle_df_frame['i'].max()
    distFunc = particle_df_frame['v_i']
    hist = distFunc.hist(ax=ax, bins=50) 
    ax.set_xlim([v_min, v_max])
    # ax.set_ylim([0, 1])
    ax.set_title(f'Frame {frame}, N = {N}')

animation = FuncAnimation(fig, update, frames=len(particle_df_list), interval=200, repeat=False)
"""


plt.show()