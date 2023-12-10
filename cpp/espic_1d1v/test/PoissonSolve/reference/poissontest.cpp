#ifndef PSTEST_HPP
#define PSTEST_HPP

#define _USE_MATH_DEFINES

#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include <cmath>
#include <fstream>
#include <cstddef>

#include "../../../include/Grid.hpp"
#include "../../../include/espic1d1v.hpp"
#include "oldpoisson.hpp" // field solver from original code

// slight refactor of field solve from original code: Vector -> std::vector<double>
size_t FieldSolveMatrix(const Eigen::SparseMatrix<double> &A, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig, const std::vector<double> &rho, std::vector<double> &phi, std::vector<double> &E_grid, const double dx);

// Driver code for testing Poisson solve
int main(int argc, char *argv[]){
    size_t Nx = std::stoi(argv[1]), routineFlag;
    double x_min = -M_PI, x_max = M_PI;
    Grid1d1v testGrid(Nx, x_min, x_max), oldGrid(Nx, x_min, x_max); // oldGrid is for *slight* refactor of original version
    double dx = testGrid.getDX();

    // This is EXACTLY how the original version solved the field
    Vector phi(Nx), rho(Nx), E_grid(Nx);

    // Initialize solver
    Eigen::SparseMatrix<double> A(Nx-1,Nx-1); // Periodic Boundary Conditions so don't need last line 
    A.reserve(Eigen::VectorXi::Constant(Nx-1,3)); // Poisson's equation so triangular
    routineFlag = BuildSparseLapl(A,dx);
    Eigen::VectorXd rhoEig(A.rows()), phiEig(A.rows());

    /* Testing should be refactored to have a wrapper around it */
    // Case 1: rho(x) = sin(x)
    std::ofstream sinLog("./data/sin.csv");
    sinLog << "x_j," << "rho_x," << "phi_x," << "phi_x_old," << "phi_original," << "E_j" << std::endl;
    for (size_t j = 0; j < Nx; j++){
        testGrid.RhoX(j) = sin(testGrid.Xgrid(j));
        oldGrid.RhoX(j) = sin(oldGrid.Xgrid(j));
        rho(j) = sin(testGrid.Xgrid(j));
    }
    routineFlag = FieldSolveMatrix(A, testGrid, rhoEig, phiEig, dx, Nx); // refactored version
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, oldGrid.getRhoX(), oldGrid.getPhiX(), oldGrid.getEX(), dx); // original version 
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, dx, rho, phi, E_grid);
    /* Stream output into csv */
    for (size_t j = 0; j < Nx; j++){
        sinLog << testGrid.Xgrid(j) << "," << testGrid.RhoX(j) << "," << testGrid.PhiX(j) << "," << oldGrid.PhiX(j) << "," << phi(j) << "," << testGrid.EX(j) << std::endl;
    }

    // Case 2: rho(x) = cos(x)
    std::ofstream cosLog("./data/cos.csv");
    cosLog << "x_j," << "rho_x," << "phi_x," << "phi_x_old," << "phi_original," << "E_j" << std::endl;
    for (size_t j = 0; j < Nx; j++){
        testGrid.RhoX(j) = cos(testGrid.Xgrid(j));
        oldGrid.RhoX(j) = cos(oldGrid.Xgrid(j));
        rho(j) = cos(testGrid.Xgrid(j));
    }
    routineFlag = FieldSolveMatrix(A, testGrid, rhoEig, phiEig, dx, Nx); // refactored version
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, oldGrid.getRhoX(), oldGrid.getPhiX(), oldGrid.getEX(), dx); // original version 
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, dx, rho, phi, E_grid);
    /* Stream output into csv */
    for (size_t j = 0; j < Nx; j++){
        cosLog << testGrid.Xgrid(j) << "," << testGrid.RhoX(j) << "," << testGrid.PhiX(j) << "," << oldGrid.PhiX(j) << "," << phi(j) << "," << testGrid.EX(j) << std::endl;
    }

    // Case 3: rho(x) = x^2
    std::ofstream paraLog("./data/parabola.csv");
    paraLog << "x_j," << "rho_x," << "phi_x," << "phi_x_old," << "phi_original," << "E_j" << std::endl;
    for (size_t j = 0; j < Nx; j++){
        testGrid.RhoX(j) = pow(testGrid.Xgrid(j),2);
        oldGrid.RhoX(j) = pow(oldGrid.Xgrid(j),2);
        rho(j) = pow(testGrid.Xgrid(j),2);
    }
    routineFlag = FieldSolveMatrix(A, testGrid, rhoEig, phiEig, dx, Nx);
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, oldGrid.getRhoX(), oldGrid.getPhiX(), oldGrid.getEX(), dx); // original version 
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, dx, rho, phi, E_grid);
    /* Stream output into csv */
    for (size_t j = 0; j < Nx; j++){
        paraLog << testGrid.Xgrid(j) << "," << testGrid.RhoX(j) << "," << testGrid.PhiX(j) << "," << oldGrid.PhiX(j) << "," << phi(j) << "," << testGrid.EX(j) << std::endl;
    }

    // Case 4: rho(x) = e^x
    std::ofstream expLog("./data/exponential.csv");
    expLog << "x_j," << "rho_x," << "phi_x," << "phi_x_old," << "phi_original," << "E_j" << std::endl;
    for (size_t j = 0; j < Nx; j++){
        testGrid.RhoX(j) = exp(testGrid.Xgrid(j));
        oldGrid.RhoX(j) = exp(oldGrid.Xgrid(j));
        rho(j) = exp(testGrid.Xgrid(j));
    }   
    routineFlag = FieldSolveMatrix(A, testGrid, rhoEig, phiEig, dx, Nx);
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, oldGrid.getRhoX(), oldGrid.getPhiX(), oldGrid.getEX(), dx); // original version 
    routineFlag = FieldSolveMatrix(A, rhoEig, phiEig, dx, rho, phi, E_grid);
    /* Stream output into csv */
    for (size_t j = 0; j < Nx; j++){
        expLog << testGrid.Xgrid(j) << "," << testGrid.RhoX(j) << "," << testGrid.PhiX(j) << "," << oldGrid.PhiX(j) << "," << phi(j) << "," << testGrid.EX(j) << std::endl;
    }

    return 0;
}
#endif

// Original version with slight refactor
size_t FieldSolveMatrix(const Eigen::SparseMatrix<double> &A, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig, const std::vector<double> &rho, std::vector<double> &phi, std::vector<double> &E_grid, const double dx){
    size_t status = 0;
    // Initialize VectorXd's
    for (size_t ij = 0; ij < rhoEig.rows(); ij++){
        rhoEig[ij] = rho[ij]; // A*phi = -rho
        phiEig[ij] = 0.0; // just initialize phi to 0 for simplicity
    }  

    // Solve A*phi = -rho for phi
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    phiEig = solver.solve(rhoEig);
    double residual = (A * phiEig - rhoEig).norm();
    // cout << "LU residual is " << residual << endl;

    // Copy out
    for (size_t ij = 0; ij < phiEig.rows(); ij++){
        phi[ij] = phiEig[ij];  
    }
    phi[phi.size() - 1] = phi[0]; // PBC

    // Calculate grid electric field 
    for (size_t ij = 0; ij < E_grid.size(); ij++){
        if (ij == 0) // first-order central difference
            E_grid[ij] = -(phi[1] - phi[phi.size()-2]) / (2.0 * dx);
        else if (ij == (E_grid.size() - 1))
            E_grid[ij] = E_grid[0]; // E_{Nx-1} = E_{0} is second PBC
            // E_grid(ij) = (phi(phi.num_rows() - 2) - phi(0)) / (2.0 * dx);
        else E_grid[ij] = -(phi[ij + 1] - phi[ij - 1]) / (2.0 * dx); 
    }
    return status;
}