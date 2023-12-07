# Overview
This version of the high-performance 1D1V C++ PIC code uses [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) to perform a sparse solve of the electrostatic field.

# Validation
Here are the cases to run to validate the code:
(1) Electron motion
- Show the correct oscillatory motion of a pair of electrons in phasespace

(2) Langmuir Waves
- Show the correct behavior of a population of electrons

(3) Twostream instability

# Directory Structure
benchmarking/

build/

data/

debug/

include/

src/

test/

# Run
**This information is outdated.** Updated information pending.

The C++ library `Eigen` is required to run this code. It was developed on v3.4.0 so this is preferred, but hopefully other versions work. 
The library must be provided to the compiler at compile-time by running the command below,

g++ -o main main.cpp -I path/to/eigen  
./main N Nx Nt W vflag

- N: number of particles
- Nx: number of grid points 
- Nt: number of timesteps
- W: weighting order (0 or 1)
- vflag: assignment flag for velocity initialization with N = 2 particles (0 or 1)

