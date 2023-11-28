# Summary
This is the C++ version of the code. 

# Directory Structure
espic_1d1v/
- Original, mixed-Eigen/STL version of the code.

cpp-petsc/
- Parallel, PETSc version of the code. 

# Current Tasks
espic_1d1v/
- Implement data collection.
- Run `parsePerfData.py` on output of kernel benchmarking on Inspiron 7591.  
- **LONG-TERM** Completely eliminate the heavily-templated `Eigen` via implementing sparse matrix solve by hand. 

cpp-petsc/
- Learn how to solve Poisson's equation on a structured grid with PETSc.

