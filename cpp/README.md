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
cpp-eigen/
- Refactor original code to work with `Eigen` containers.

cpp-petsc/
- Need to learn PETSc.

espic_1d1v/
- Clean up, and refactor code.
- Build with compiler optimizations, and benchmark.
- Completely eliminate the slow, heavily-templated `Eigen` via implementing sparse matrix solve by hand. 