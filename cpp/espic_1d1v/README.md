# Run Instructions:
# Overview
This version of the high-performance 1D1V C++ PIC code uses [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) to perform a sparse solve of the electrostatic field, as mentioned in the README of the root directory for this repository. 

As a consequence of this, an additional amount of O(Nx) work must be done every timestep in moving the data from the `std::vector` STL container into the `VectorXd` `Eigen` container. This somewhat defeats the reason to use a sparse solve in the first place, namely to achieve performance at fine grids (large Nx); however, there are several other places in the program where O(Nx) work must be performed, so ultimately this is not a big deal. That being said, care must be taken when working with `Eigen` to not destroy the performance with a bad allocation strategy.   

# Run
The C++ library `Eigen` is required to run this code. It was developed on v3.4.0 so this is preferred, but hopefully other versions work. 
The library must be provided to the compiler at compile-time by running the command below,

g++ -o main main.cpp -I path/to/eigen  
./main N Nx Nt W vflag

- N: number of particles
- Nx: number of grid points 
- Nt: number of timesteps
- W: weighting order (0 or 1)
- vflag: assignment flag for velocity initialization with N = 2 particles (0 or 1)

