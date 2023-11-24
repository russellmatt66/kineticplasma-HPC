# pic1d1v-HPC
This repository houses a high-performance Particle-In-Cell (PIC) code that simulates a plasma from a kinetic perspective. This is in contrast to a fluid point of view where the plasma is viewed as a continuum that is governed explicitly by a set of conservation laws. For an explanation of the PIC algorithm, see the texts, *Plasma Physics via Computer Simulation* or *Computer Simulation Using Particles*. The file located in this directory, "The-PIC-Algorithm.pdf", describes the code that is implemented here. 

# Directory Structure
cpp/
- Work is currently in progress in `cpp/espic_1d1v/` to refactor the original version of the code which utilizes the `Eigen` library to perform a sparse solve of the electrostatic field. 

c/
- Development on the C version of the code has not yet begun.

julia/
- Development of the Julia version of the code has not yet begun.

python/
- Development of the Python version of the code has not yet begun.

fortran/
- Development of the Fortran version of the code has not yet begun.

rust/ 
- Development of the Rust version of the code has not yet begun.

# A Competition Among Programming Languages 
In addition to the original code implemented in C++, versions of the code implemented in Fortran, Python, Rust, Julia, and C are being developed in the respective folders for the purpose of comparing the performance of the different languages. 