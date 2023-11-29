# Summary
This is the C++ version of the code. 

# Directory Structure
espic_1d1v/
- Original, mixed-Eigen/STL version of the code.

cpp-petsc/
- Parallel, PETSc version of the code. 
- Development has not yet begun on this version of the code.

twostream/
- Application code that uses the PIC kernel developed in `espic_1d1v/` to simulate the two-stream instability.

# Current Tasks
espic_1d1v/
- Implement data collection from simulation.
- Continue to analyze output data from kernel benchmarking.
- Write up relevant parts of master report

cpp-petsc/
- Learn how to solve Poisson's equation on a structured grid with PETSc.

twostream/
- Begin developing this project.
