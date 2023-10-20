# Run Instructions:
The C++ library 'Eigen' is required to run this code. It was developed on v3.4.0 so this is preferred, but hopefully other versions work. 
The library must be provided to the compiler at compile-time by running the command below,

g++ -o main main.cpp -I path/to/eigen  
./main N Nx Nt W vflag

- N: number of particles
- Nx: number of grid points 
- Nt: number of timesteps
- W: weighting order (0 or 1)
- vflag: assignment flag for velocity initialization with N = 2 particles (0 or 1)

# Overview
This is the version of the high-performance 1D1V PIC code that uses Eigen to perform a sparse solve of the electrostatic field, as mentioned in the README of the root directory for this repository. To read a description of the PIC algorithm, look inside the file "SP23\_AA545\_ComputerProject1_2.pdf". 

This code is largely obsolete, and should not be considered anything more than a proof-of-concept. It possesses no build system, takes command-line input instead of a .inp file, uses .csv files to collect the data, and a mountain of Python spaghetti to process it. For a modern version, that is currently under development, look inside the `cpp-petsc/` directory. 
