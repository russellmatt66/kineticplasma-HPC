#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include<cmath>
#include<Eigen/Sparse>
#include<random>

#include "test.hpp"
#include "../espic.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    size_t N = stoul(argv[1]); // number of particles
    size_t Nx = stoul(argv[2]); // number of gridpoints (# grid cells = Nx - 1)
    size_t Nt = stoul(argv[3]); // number of timesteps
    size_t W = stoul(argv[4]); // particle/force-weighting order (0th or 1st)
    size_t vflag = stoul(argv[5]); // only use with N = 2, flag for whether to

    Vector x_grid = Vector(Nx), rho = Vector(Nx), phi = Vector(Nx), E_grid = Vector(Nx); 
    double x_min = -M_PI, x_max = M_PI;
    double dx = ( x_max - x_min ) / (x_grid.num_rows() - 1);
    for (size_t ij = 0; ij < x_grid.num_rows(); ij++){
        x_grid(ij) = x_min + ij * dx;
        rho(ij) = 0.0;
    }
   
    Vector particle_pos(N);
    particle_pos(0) = -M_PI / 4.0;
    particle_pos(1) = M_PI / 4.0;

    ofstream testgrid,testenergy,testparticles;
    testgrid.open("test_grid.csv");
    testenergy.open("test_energy.csv");
    testparticles.open("testparticles.csv");
    
    // initialize
    size_t routineFlag = ParticleWeightTest(W, particle_pos, x_grid, rho);

    // loop
    for (size_t it = 1; it < N; it++){
        size_t routineFlag = ParticleWeightTest(W,particle_pos,x_grid,rho);
        size_t routineFlag = FieldSolveTest(A,dx,rho,phi,E_grid);
        size_t routineFlag = ForceWeightTest(W,E_grid,,x_grid,,particle_x,particle_E);
        size_t routineFlag = ParticlePushTest();
    }
    return 0;
}
