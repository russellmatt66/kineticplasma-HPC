import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy

def parseParticles(particle_df: pd.DataFrame) -> tuple: 
    """
    @brief: 
    Works through the DataFrame containing the particle output data and places the data associated with each particle, i, into
    a list, sorted in ascending order w.r.t the particle number, i.
    @input: 
    particle_df - the DataFrame containing the raw data from the .csv file.
    @output:
    list_particle_df - a list of DataFrames, each containing the data for particle i, sorted in ascending order. 
    N - number of particles
    color_list - list of colors to use for plotting different particles
    """
    list_particle_df = []
    N = particle_df['i'].max() + 1 # number of particles - starts at 0
    Nt = particle_df['n'].max() + 1 # number of timesteps - starts at 0
    for ip in range(N):
        temp_df = particle_df[particle_df['i'] == ip]
        list_particle_df.append(temp_df)
    # for N particles, need N colors
    cmap = plt.cm.get_cmap('plasma', N)
    colors = [cmap(i) for i in range(N)]
    hex_strings = [mcolors.to_hex(color) for color in colors] # '#rrggbb'
    hex_colors = [hex_str[:].upper() for hex_str in hex_strings] # '#RRGGBB' 
    # print(hex_colors)

    return list_particle_df, N, hex_colors, Nt

def plotPhaseSpace(list_particle_df: list, N: numpy.int64, hex_colors: list, phaseSpaceAx: plt.Axes) -> numpy.int64:
    """
    @brief: plots trajectory of the particles in phase-space
    @input:
    (1) list_particle_df - the list of DataFrames containing the particle trajectories
    (2) N - the number of particles
    (3) hex_colors - a list of hexadecimal codes corresponding to the color of the corresponding particle 
    (4) phaseSpaceAx - the Axes object the trajectories are going to be plotted on
    @output: 
    status - just a status flag
    """
    status = numpy.int64(1)
    # Code goes here - prototype in parse_and_plot.py
    for ip in range(N):
        # print(list_particle_df[ip])
        # print(type(list_particle_df[ip]))
        particles_df_ip = list_particle_df[ip] 
        # particles_df_ip.plot(kind='scatter', x='position', y = 'velocity', ax=phaseSpaceAx, color=hex_colors[ip],label='Particle {}'.format(ip)) 
        particles_df_ip.plot(kind='scatter', x='position', y = 'velocity', ax=phaseSpaceAx, color=hex_colors[ip]) 
    return status

def zeroCrossings(energy_history_df: pd.DataFrame) -> int: # Get number of cycles
    average_E = energy_history_df['Total Energy'].mean()
    energy_history_df['shifted energy'] = energy_history_df['Total Energy'] - average_E
    sign_changes = np.sign(energy_history_df['shifted energy'].diff().fillna(0)).diff()
    zerocrossings = np.where(sign_changes > 0)[0]
    return len(zerocrossings)

# def makeGridMovie():