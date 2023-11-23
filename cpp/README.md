# Summary
This is the C++ version of the code. 

# Directory Structure
cpp-eigen/
- Sequential, Eigen version of the code.

cpp-petsc/
- Parallel, PETSc version of the code. 

espic_1d1v/
- Original, mixed-Eigen/STL version of the code. Currently being refactored.

# Current Tasks
espic_1d1v/
- Implement integration of charge across grid.
- Implement initial conditions for species.
- Implement data collection.
- Build with compiler optimizations, and benchmark with hardware performance counters.
- **LONG-TERM** Completely eliminate the slow, heavily-templated `Eigen` via implementing sparse matrix solve by hand. 

cpp-eigen/
- Refactor original code to be based on `Eigen` containers.

cpp-petsc/
- Learn how to solve Poisson's equation on a structured grid with PETSc.

