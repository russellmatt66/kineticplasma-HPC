# pic1d1v-HPC
This repository houses a high-performance Particle-In-Cell (PIC) code that simulates a plasma from a kinetic perspective. This is in contrast to a fluid point of view where the plasma is viewed as a continuum that is governed explicitly by a set of conservation laws. 

"1d1v" refers to the size of the phasespace that the simulation particles are sampled from: one spatial dimension (1d), and one velocity dimension (1v). This is the minimum relevant configuration for a PIC code. 

Performant code is achieved by implementing a binary search for finding the particles, and a sparse matrix solver for computing the field. 

A working version of this code can be found in the uw_MSAA-classwork repository. It is currently in the process of being ported over to here. Additionally, versions of the code implemented in Fortran, Python, and C are being developed to compare the performance of the different languages.
