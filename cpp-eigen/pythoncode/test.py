import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import helper

particles_df_test = pd.read_csv('particles_phasespacetest.csv')
grid_df_test = pd.read_csv('grid_datatest.csv')
energy_df_test = pd.read_csv('energy_historytest.csv')

""" Phasespace plotting"""
psfig_test, psax_test = plt.subplots() 
particles_df_test_list, N_test, particlecolors_Ntest, Nt_test = helper.parseParticles(particles_df_test)
vprime = particles_df_test.loc[particles_df_test['t']==0, 'velocity'].iloc[0]
routineFlag = helper.plotPhaseSpace(particles_df_test_list, N_test, particlecolors_Ntest, psax_test)
psfig_test.suptitle('Phase-space trajectories for $N = %s$, $N_{t} = %s$' %(N_test,Nt_test))
plt.legend()

""" Energy History Plotting """
Efig_test, Eax_test = plt.subplots()
energy_df_test.plot(kind='scatter', x = 't' , y = 'Kinetic Energy', ax=Eax_test, color='blue', label = 'Kinetic Energy')
energy_df_test.plot(kind='scatter', x = 't', y = 'Electric Energy', ax=Eax_test, color='red', label = 'Electric Energy')
energy_df_test.plot(kind='scatter', x = 't', y = 'Total Energy', ax=Eax_test, color='green', label = 'Total Energy')
Efig_test.suptitle('Energy history for $N = %s$, $N_{t} = %s$' %(N_test,Nt_test))
plt.legend()

plt.show()