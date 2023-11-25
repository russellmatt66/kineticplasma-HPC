# Summary
This folder builds an executable which can be used to benchmark the kernel of the PIC code. Additionally, a bash script 'perf_data.sh', and a python script 'CollectPerfData.py', are used to benchmark the code by counting hardware floating-point events.

# Dell Inspiron 7591 2n1
N_max = 2^28 = 268435456
Nx_max = 2^13 = 8192

# Build
`Eigen` must be installed, and CMake must be pointed to the folder where the `FindEigen3.cmake` file is. 

On my machine, this is accomplished via,

    cmake ../src -DCMAKE_MODULE_PATH=$HOME/software/eigen-3.4.0/cmake

Generally, it will be of the form,

    cmake ../src -DCMAKE_MODULE_PATH=$HOME/path/to/eigen-3.4.0/cmake

# Run
Configure the input file `bench.inp` with the appropriate values, and then execute the binary using the appropriate hardware code to count the floating point events:
    
    perf stat -e r5301c7 ./espic1d1v-bench

is the command on Skylake architecture.

- N: number of particles to simulate
- Nx: number of gridpoints
- Nt: number of timesteps
- dtcoeff: The fraction of a whole plasma period that one timestep is
- W: Weighting order (related to particle shape)
- vprime: Amplitude of velocity perturbation, or related to FWHM of random pop
- particleICs: How the particle population should be initialized. Either 'random' or 'uniform'

# Analyze
