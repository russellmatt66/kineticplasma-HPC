# pic1d1v-HPC
This repository houses a high-performance Particle-In-Cell (PIC) code that simulates a plasma from a kinetic perspective. This is in contrast to a fluid point of view where the plasma is viewed as a continuum that is being governed by magnetohydrodynamic (MHD) equations. Performant code is achieved by implementing a binary search for finding the particles, and a sparse matrix solver for computing the field. 

"1d1v" refers to the size of the phasespace that the simulation particles are sampled from: one spatial dimension (1d), and one velocity dimension (1v). This is the minimum relevant configuration for a PIC code. 

A working version of this code can be found in the uw_MSAA-classwork repository. It is currently in the process of being ported over to here. 
