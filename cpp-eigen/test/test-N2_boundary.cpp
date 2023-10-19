#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include<cmath>
#include<Eigen/Sparse>
#include<random>

#include "test.hpp"
#include "../espic.hpp"

using namespace std;

// Testing what happens when particle is near boundary
int main(int argc, char* argv[])
{
    size_t Nt = stoul(argv[1]);

    size_t Nx = 33;
    Vector x_grid = Vector(Nx), rho = Vector(Nx), phi = Vector(Nx), E_grid = Vector(Nx);
    double x_min = -M_PI, x_max = M_PI, dx = (x_max - x_min) / (Nx - 1), Q_net;
    
    // Initialize grid
    for(size_t ij = 0; ij < Nx; ij++){
        x_grid(ij) = x_min + ij * dx;
        rho(ij) = 0.0;
        phi(ij) = 0.0;
        E_grid(ij) = 0.0;
    }

    // Initialize particles
    size_t N = 2, W = 1, routineFlag;
    Vector particle_x = Vector(N), particle_v = Vector(N), particle_E = Vector(N), particle_indices = Vector(N);
    particle_x(0) = -dx/10.0;
    particle_x(1) = M_PI - dx / 10.0;
    particle_v(0) = 0.0;
    particle_v(1) = 0.0

    // Initialize time
    double omega_p = sqrt(N/(x_max - x_min)); // sqrt(N/L)
    double tau_p = (2.0 * M_PI) / omega_p;
    double t = 0.0; // long doubles are a nightmare
    double dt = 0.01 * tau_p; // should this be tau_p or omega_p?


    
    return 0;
}
