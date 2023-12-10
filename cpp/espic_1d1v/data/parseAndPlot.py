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
energy_path = 'Energy'

# Creates a list of all relevant files based on the argument
grid_files = glob.glob(grid_path + "/*.csv")
particle_files = glob.glob(particle_path + "/*.csv")

# Little slower than using generator expression, i.e., () instead of []
grid_df_list = [pd.read_csv(datafile) for datafile in grid_files]
particle_df_list = [pd.read_csv(datafile) for datafile in particle_files]

# print(particle_df_list.type())

grid_df = pd.concat(grid_df_list, ignore_index=True)
particle_df = pd.concat(particle_df_list, ignore_index=True)
energy_df = pd.read_csv(energy_path + "/EnergyHistory.csv")

# Looks good here
# print(grid_df)
# print(particle_df)

"""
Parse DataFrames
"""


"""
Plot output
"""
# Create plot of energy history
energyFig, energyAx = plt.subplots()

energyAx.plot(energy_df['n'],energy_df['KE'])
energyAx.plot(energy_df['n'],energy_df['PE'])
energyAx.plot(energy_df['n'],energy_df['E'])

# Create movie of grid field and potential


# Create plot of phase-space
phaseFig, phaseAx = plt.subplots()

ParticleX = particle_df['x_i']
ParticleVx = particle_df['v_i']

phaseAx.scatter(ParticleX, ParticleVx)

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