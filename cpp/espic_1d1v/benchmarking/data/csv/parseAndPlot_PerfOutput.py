import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Read CSV file
perfoutput = pd.read_csv("PerfOutput_Kernel.csv")

NumArith = perfoutput.FP_ARITH.str.replace(',','') # perf stat -e provides r5301c7 number with commas

# Performance = pd.to_numeric(NumArith) / perfoutput.WALLTIME_MEAN # Related to average performance

AvgPerformance = (pd.to_numeric(NumArith) / 25) * perfoutput.INVWALLTIME_SUM # 25 is a magic number corresponding to number of trials, care should be taken in the future to store this number in csv

# Find where max and minimum performance is achieved
maxPerf = AvgPerformance.max() / 1e9
minPerf = AvgPerformance.min() / 1e9

maxInd = AvgPerformance.argmax()
minInd = AvgPerformance.argmin()

N_max = perfoutput.N[maxInd]
Nx_max = perfoutput.Nx[maxInd]
N_min = perfoutput.N[minInd]
Nx_min = perfoutput.Nx[minInd]

# Find where performance is above 1 Gflop, 2 Gflop
sorted_perfoutput = perfoutput.sort_values(by=['N', 'Nx'])
sorted_NumArith = sorted_perfoutput.FP_ARITH.str.replace(',','')
sorted_AvgPerformance = (pd.to_numeric(sorted_NumArith) / 25) * sorted_perfoutput.INVWALLTIME_SUM

gt1Gflop = sorted_AvgPerformance[(sorted_AvgPerformance > 1e9) & (sorted_AvgPerformance < 2e9)]
gt1Gflop = gt1Gflop.index
gt2Gflop = sorted_AvgPerformance[sorted_AvgPerformance > 2e9].index

print(gt1Gflop)
print(gt2Gflop)

gt1Gflop_N = sorted_perfoutput.N[gt1Gflop]
gt1Gflop_Nx = sorted_perfoutput.Nx[gt1Gflop]

gt2Gflop_N = sorted_perfoutput.N[gt2Gflop]
gt2Gflop_Nx = sorted_perfoutput.Nx[gt2Gflop]

print(gt1Gflop_N)
print(gt1Gflop_Nx)
print(gt2Gflop_N)
print(gt2Gflop_Nx)

with open('maxMin.txt', 'w', newline = '') as fp:
    fp.write("Maximum performance achieved is " + str(maxPerf) + " Gflops, with N=" + str(N_max) + ",and Nx=" + str(Nx_max) + "\n")
    fp.write("Minimum performance achieved is " + str(minPerf) + " Gflops, with N=" + str(N_min) + ",and Nx=" + str(Nx_min) + "\n")
    fp.write("\n")
    fp.write("Performance greater than 2 Gflops is achieved for the following problem sizes: \n")
    for index in gt2Gflop:
        fp.write("(N,Nx) = " + "(" + str(perfoutput.N[index]) + "," + str(perfoutput.Nx[index]) + ")\n")
    fp.write("\n")
    fp.write("Performance between 1 Gflops and 2 Gflops is achieved for the following problem sizes: \n")
    for index in gt1Gflop:  
        fp.write("(N,Nx) = " + "(" + str(perfoutput.N[index]) + "," + str(perfoutput.Nx[index]) + ")\n") 


# Plot data "landscapes"
fig = plt.figure()
ax = fig.add_subplot(211, projection='3d')

# Plot landscape of the Mean Walltime
ax.plot_trisurf(np.log2(perfoutput.N), np.log10(perfoutput.Nx), perfoutput.WALLTIME_MEAN, cmap='viridis')

# Set labels
ax.set_xlabel('$log_{2}$(N)')
ax.set_ylabel('$log_{10}$(Nx)')
ax.set_zlabel('$\\bar{\\tau}_{wall}$ (s)')
ax.set_title('Average Walltime of 1D1V ES-PIC Simulation Kernel')


# # Plot landscape of the Kernel "Performance" 
# ax = fig.add_subplot(222, projection='3d')
# ax.plot_trisurf(np.log2(perfoutput.N), np.log10(perfoutput.Nx), Performance / 1e9, cmap='viridis')

# # Set labels
# ax.set_xlabel('$log_{2}$(N)')
# ax.set_ylabel('$log_{10}$(Nx)')
# ax.set_zlabel('Performance (Gflops)')
# ax.set_title('Performance of 1D1V ES-PIC Simulation Kernel')

# fig.suptitle("Benchmarking a Serial 1D1V ES-PIC Simulation Kernel on a Dell Inspiron 7591 2n1")

# Plot landscape of the Kernel Mean Performance
ax = fig.add_subplot(212, projection='3d')
ax.plot_trisurf(np.log2(perfoutput.N), np.log10(perfoutput.Nx), AvgPerformance / 1e9, cmap='viridis')

# Set labels
ax.set_xlabel('$log_{2}$(N)')
ax.set_ylabel('$log_{10}$(Nx)')
ax.set_zlabel('Average Performance (Gflops)')
ax.set_title('Average performance of 1D1V ES-PIC Simulation Kernel')

fig.suptitle("Benchmarking a Serial 1D1V ES-PIC Simulation Kernel on a Dell Inspiron 7591 2n1")

plt.show()