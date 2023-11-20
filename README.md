# pic1d1v-HPC
This repository houses a high-performance Particle-In-Cell (PIC) code that simulates a plasma from a kinetic perspective. This is in contrast to a fluid point of view where the plasma is viewed as a continuum that is governed explicitly by a set of conservation laws. For an explanation of the PIC algorithm, see the file located in this directory, "The-PIC-Algorithm.pdf". "1D1V" refers to the size of the phasespace that the simulation particles are sampled from: one spatial dimension (1d), and one velocity dimension (1v). This is the minimum relevant configuration for a PIC code. 

# High-Performance: More than a Buzzword
Performant code is achieved by implementing a binary search for finding the particles, and a sparse solver for computing the fields, which allows problem sizes of >1,000,000 particles to be simulated with a sequential performance of ~200 Mflops without any compiler optimizations. The original version of the code, located in `cpp-eigen/`, uses the Eigen library to implement this, however, in this situation the performance does not scale with increasing grid fineness because of the need to perform O(Nx) work every timestep moving data between `std::vector` and the Eigen container needed by the sparse solver. 

In order to address this, the code was refactored to base all data containers on the Eigen library. This addressed the issue of the performance degradation at fine grids by degrading the performance across the board to ~10 Mflops. Profiling with `gprof` showed that this was due to the introduction of a large amount of overhead as a result of using the Eigen library. The branch where this work was done is not located in this repository.

To achieve performant code at fine grids, and for large numbers of particles, a final refactor based on utilizing the PETSC library is occurring in the `cpp-petsc/` directory.

# A Competition Among Programming Languages 
Additionally, versions of the code implemented in Fortran, Python, Rust, Julia, and C are being developed in the respective folders for the purpose of comparing the performance of the different languages. The expectation is that C will be the highest-performing language, but the author's varying degree of skill with the different languages leaves things up in the air.
