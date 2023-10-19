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
        rho(ij) = sin(x_grid(ij));
        phi(ij) = 0.0;
        E_grid(ij) = 0.0;
    }

    ofstream grid_file;
    grid_file.open("fieldsolvetest_grid_ansoln.csv"); 
    grid_file << "j" << "," << "x_grid(j)" << "," << "phi" << "," << "rho" << "," << "E_grid" << endl;

    // Initialize Laplacian Matrix Operator
    Eigen::SparseMatrix<double> Lapl(Nx,Nx);
    Lapl.reserve(Eigen::VectorXi::Constant(Nx,3));
    routineFlag = BuildSparseLapl(Lapl,dx);

    // Initialize First Derivative Matrix Operator
    // Eigen::SparseMatrix<double> FD(Nx,Nx);
    // FD.reserve(Eigen::VectorXi::Constant(Nx,3));
    // routineFlag = BuildSparseFD(FD,dx);

    routineFlag = FieldSolveMatrix(Lapl,dx,rho,phi,E_grid);

    for (size_t ij = 0; ij < Nx; ij++){
        grid_file << ij << "," << x_grid(ij) << "," << phi(ij) << "," << rho(ij) << "," << E_grid(ij) << endl;
    }
    
    grid_file.close();
    return 0;
}