#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include<cmath>
#include<Eigen/Sparse>
#include<random>

#include "test.hpp"
#include "../espic.hpp"

using namespace std;

int main()
{   
    size_t N = 64;
    Vector particle_x = Vector(N), particle_v = Vector(N);
    double x_min = -M_PI, x_max = M_PI;
    if (N == 64){ // sinusoidal perturbation
        double k = 1.0;
        double dx_particles = (x_max - x_min) / (N - 1);
        cout << "Particle spacing " << dx_particles << endl;
        for (size_t ii = 0; ii < particle_x.num_rows(); ii++){
            particle_x(ii) = x_min + double(ii) * dx_particles;
            particle_v(ii) = sin(k*(particle_x(ii)));
            cout << "Particle " << ii << " at " << particle_x(ii) << endl;
        }
        
    }
    

    // Initialize grid 
    // size_t Nx = 33;
    // Vector x_grid = Vector(Nx), rho = Vector(Nx), phi = Vector(Nx), E_grid = Vector(Nx);
    // double x_min = -M_PI, x_max = M_PI;
    // double dx = ( x_max - x_min ) / (x_grid.num_rows() - 1);
    // for (size_t ij = 0; ij < x_grid.num_rows(); ij++)
    //     x_grid(ij) = x_min + ij * dx; // uniform grid
    
    // Eigen::SparseMatrix<double> A(Nx,Nx);
    // A.reserve(Eigen::VectorXi::Constant(Nx,3));
    // size_t routineFlag = BuildSparseMatrix(A,dx);
    
    // // Test - equidistant particles = no electric field force
    // size_t N = 2;
    // size_t W = 1;
    // double Q_net;
    // Vector particle_x = Vector(N), particle_v = Vector(N), particle_E = Vector(N);
    // particle_x(0) = -M_PI / 4.0;
    // particle_x(1) = M_PI / 4.0;
    // particle_v(0) = 0.0;
    // particle_v(1) = 0.0;

    // // Charge-weight
    // routineFlag = ParticleWeightTest(W,particle_x,x_grid,rho,Q_net);
    // cout << "Charge weighted " << endl;
    // // Field-solve
    // routineFlag = FieldSolve(A,dx,rho,phi,E_grid);
    // cout << "Field solved " << endl;
    // // Force-weight
    // routineFlag = ForceWeight(W,E_grid,x_grid,particle_x,particle_E);

    // cout << "Electric force on the particles " << endl;
    // for (size_t ii = 0; ii < particle_E.num_rows(); ii++){
    //     cout << particle_E(ii) << ",";
    // }
    // cout << endl;

    // Test - no particles = no field
    // for (size_t ij = 0; ij < rho.num_rows(); ij++){
    //     rho(ij) = 0.0;
    // }

    // routineFlag = FieldSolve(A,dx,rho,phi,E_grid);

    // cout << "Electric field " << endl;
    // for (size_t ij = 0; ij < E_grid.num_rows(); ij++)
    //     cout << E_grid(ij);
    // cout << endl;

    // Test - Poisson Solver Stencil
    // size_t Nx = 33;
    // Eigen::SparseMatrix<double> A(Nx,Nx);

    // for (size_t iA = 0; iA < Nx; iA++){
    //     if (iA > 0) A.insert(iA,iA-1) = 1.0;
    //     A.insert(iA,iA) = -2.0;
    //     if (iA < Nx - 1) A.insert(iA,iA+1) = 1.0; 
    // }
    // // Periodic Boundary Conditions
    // // phi(x) = phi(x+L) -> (j = 0): phi_{0} - phi_{Nx-1} = 0 
    // A.coeffRef(0,0) = 1.0;
    // A.coeffRef(0,1) = 0.0;
    // A.coeffRef(0,Nx-1) = -1.0; 
    // // A.coeffRef(0,Nx-1) = 0.0;
    
    // // dphi/dx|x=x_0 = dphi/dx|x=x_L
    // // -> (j = Nx-1): phi_{1} - phi_{0} - phi_{Nx-1} + phi_{Nx-2} = 0
    // // A.coeffRef(Nx-1,0) = -1.0;
    // A.coeffRef(Nx-1,0) = 1.0;
    // // A.coeffRef(Nx-1,1) = 1.0;
    // A.coeffRef(Nx-1,1) = 0.0;
    // A.coeffRef(Nx-1,Nx-1) = 0.0;
    // // A.coeffRef(Nx-1,Nx-1) = 1.0; // set gauge with phi_{Nx-1}
    // // A.coeffRef(Nx-1,Nx-2) = 1.0;
    // A.coeffRef(Nx-1,Nx-2) = 0.0;

    // cout << A << endl;

    // Test - Particle Weighting
    // cout << "Testing particle weighting" << endl;
    // size_t Nx = 33;
    // Vector x_grid = Vector(Nx), rho = Vector(Nx); 
    // double x_min = -M_PI, x_max = M_PI;
    // double dx = ( x_max - x_min ) / (x_grid.num_rows() - 1);
    // for (size_t ix = 0; ix < x_grid.num_rows(); ix++){
    //     x_grid(ix) = x_min + ix * dx; 
    //     cout << x_grid(ix) << endl;
    //     rho(ix) = 0.0;
    // }
    // size_t W = 1; // weighting-order
    // size_t N = 2; // number of particles
    // Vector particle_pos(N);
    // particle_pos(0) = -M_PI / 4.0;
    // particle_pos(1) = M_PI / 4.0;

    // size_t routineFlag = ParticleWeightTest(W,particle_pos,x_grid,rho);
    
    // Test - Binary Search
    // size_t Nx = 33;
    // size_t N = 64;
    // Vector x_grid = Vector(Nx), particle_pos = Vector(N); 
    // double x_min = -M_PI, x_max = M_PI;
    // double dx = ( x_max - x_min ) / (x_grid.num_rows() - 1);
    // for (size_t ix = 0; ix < x_grid.num_rows(); ix++){
    //     x_grid(ix) = x_min + ix * dx; 
    //     cout << x_grid(ix) << endl;
    // }
    
    // Generate test data
    // random_device rd;
    // mt19937 gen(rd());
    // uniform_real_distribution<double> dis(x_min, x_max);
    // for (size_t ii = 0; ii < particle_pos.num_rows(); ii++){
    //     particle_pos(ii) = dis(gen);
    // }

    // size_t j_found, j_raster;
    // ofstream bstfile; // binary search test data
    // bstfile.open("test/bst.csv"); 
    // bstfile << "i,x_i,j_found,j_raster,x_jfound" << endl;
    // for (size_t ii = 0; ii < particle_pos.num_rows(); ii++){
    //     cout << "Searching for particle #" << ii << endl;
    //     j_found = findParticleTest(x_grid,particle_pos(ii));
    //     j_raster = findParticleRaster(x_grid,particle_pos(ii));
    //     bstfile << ii << "," << particle_pos(ii) << "," << j_found << "," << j_raster << "," << x_grid(j_found) << endl;
    // }
    // bstfile.close();

    // // Test - Passing ofstream
    // size_t test_size = 100;
    // ofstream testFile;
    // Vector testData(test_size);
    // for (size_t it = 0; it < test_size; it++){ // initialize test data
    //     testData(it) = it; // implicit double conversion
    // }

    // size_t routineFlag = CollectDataTest(testFile,testData);

    // // Test - Sparse Matrix Creation
    // size_t Nx = 33;
    // Eigen::SparseMatrix<double> A(Nx,Nx);
    // A.reserve(Eigen::VectorXi::Constant(Nx,3));
    
    // // Tridiagonal stencil 
    // for (size_t iA = 0; iA < Nx; iA++){
    //     if (iA > 0) A.insert(iA,iA-1) = 1.0;
    //     A.insert(iA,iA) = -2.0;
    //     if (iA < Nx - 1) A.insert(iA,iA+1) = 1.0; 
    // }
    // // Periodic Boundary Conditions
    // // phi(x) = phi(x+L) -> phi_{0} - phi_{Nx-1} = 0
    // A.coeffRef(0,0) = 1.0;
    // A.coeffRef(0,Nx-1) = 0.0;
    // A.coeffRef(0,1) = 0.0;
    
    // // dphi/dx|x=x_0 = dphi/dx|x=x_L
    // // -> phi_{1} - phi_{0} - phi_{Nx-1} + phi_{Nx-2} = 0
    // A.coeffRef(Nx-1,0) = -1.0;
    // A.coeffRef(Nx-1,1) = 1.0;
    // A.coeffRef(Nx-1,Nx-1) = -1.0;
    // A.coeffRef(Nx-1,Nx-2) = 1.0;

    // // Print A
    // std::cout << A << std::endl; 
    
    return 0;
}
