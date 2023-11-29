import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Read CSV file
perfoutput = pd.read_csv("PerfOutput_Kernel.csv")

NumArith = perfoutput.FP_ARITH.str.replace(',','')

Performance = pd.to_numeric(NumArith) / perfoutput.WALLTIME_MEAN

# Plot data "landscapes"
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')

# Plot landscape of the Mean Walltime
ax.plot_trisurf(np.log2(perfoutput.N), np.log10(perfoutput.Nx), perfoutput.WALLTIME_MEAN, cmap='viridis')

# Set labels
ax.set_xlabel('log2(N)')
ax.set_ylabel('log10(Nx)')
ax.set_zlabel('Mean Walltime')
ax.set_title('Average Walltime of 1D1V ES-PIC Simulation Kernel')


# Plot landscape of the Kernel Performance 
ax = fig.add_subplot(122, projection='3d')
ax.plot_trisurf(np.log2(perfoutput.N), np.log10(perfoutput.Nx), Performance / 1e9, cmap='viridis')

# Set labels
ax.set_xlabel('log2(N)')
ax.set_ylabel('log10(Nx)')
ax.set_zlabel('Performance (Gflops)')
ax.set_title('Performance of 1D1V ES-PIC Simulation Kernel')

fig.suptitle("Benchmarking a Serial 1D1V ES-PIC Simulation Kernel on a Dell Inspiron 7591 2n1")

plt.show()