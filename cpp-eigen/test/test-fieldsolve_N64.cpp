#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<random>
#include<chrono>
#include <Eigen/Sparse>

#include "../espic.hpp"
#include "../Vector.hpp"


int main(int argc, char* argv[]){
    size_t Nx = 33, routineFlag;
    Vector x_grid = Vector(Nx), rho = Vector(Nx), phi = Vector(Nx), E_grid = Vector(Nx);
    double x_min = -M_PI, x_max = M_PI, dx = (x_max - x_min) / (Nx -1), Q_net;
    // Initialize grid
    for(size_t ij = 0; ij < Nx; ij++){
        x_grid(ij) = x_min + ij * dx;
        rho(ij) = 0.0;
        phi(ij) = 0.0;
        E_grid(ij) = 0.0;
    }
    // Initialize Laplacian Matrix Operator
    Eigen::SparseMatrix<double> Lapl(Nx,Nx);
    Lapl.reserve(Eigen::VectorXi::Constant(Nx,3));
    routineFlag = BuildSparseLapl(Lapl,dx);

    // Initialize First Derivative Matrix Operator
    // Eigen::SparseMatrix<double> FD(Nx,Nx);
    // FD.reserve(Eigen::VectorXi::Constant(Nx,3));
    // routineFlag = BuildSparseFD(FD,dx);

    // Initialize equidistant particles
    size_t N = 64, W = 1;
    Vector particle_x = Vector(N), particle_E = Vector(N), particle_indices = Vector(N);
    double dx_p = (x_max - x_min) / N, dist, threshold;
    double m = (1.0 / (N-1)) * (x_max - x_min - dx_p), b = x_min + dx_p / 2.0;
    for(size_t ii = 0; ii < N; ii++){
        particle_x(ii) = m * ii + b;
    }    

    ofstream particle_file, grid_file;
    particle_file.open("fieldsolvetest_particles.csv");
    particle_file << "i" << "," << "position" << "," << "j_found" << "," << "x_grid(j_found)" << "," << "x_grid(j_found+1)" << "," << "dist" << "," << "particle_E" << endl;
    grid_file.open("fieldsolvetest_grid.csv"); 
    grid_file << "j" << "," << "x_grid(j)" << "," << "rho" <<  "," << "phi" << "," << "E_grid" << "," << "Q_net" << endl;
    routineFlag = ParticleWeight(W, particle_x, x_grid, rho, particle_indices, Q_net);
    routineFlag = FieldSolveMatrix(Lapl, dx, rho, phi, E_grid);
    routineFlag = ForceWeight(W, particle_indices, E_grid, x_grid, particle_x, particle_E);
    for (size_t ii = 0; ii < N; ii++){
            particle_file << ii << "," << particle_x(ii) << "," << particle_indices(ii) << "," << x_grid(particle_indices(ii)) << "," << x_grid(particle_indices(ii)+1) << "," << dist << "," << particle_E(ii) << endl;
        }
    for (size_t ij = 0; ij < Nx; ij++){
        grid_file << ij << "," << x_grid(ij) << "," << rho(ij) << "," << phi(ij) << "," << E_grid(ij) << "," << Q_net << endl;
    }
    
    particle_file.close();
    grid_file.close();
    return 0;
}