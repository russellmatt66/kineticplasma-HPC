#ifndef OLD_POISSON_HPP
#define OLD_POISSON_HPP

#include <cstddef>
#include <vector>

#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include "../../include/Grid.hpp"
#include "../../include/espic1d1v.hpp"

class Vector {
public:
  Vector(size_t M) : num_rows_(M), storage_(num_rows_) {}

        double& operator()(size_t i)       { return storage_[i]; }
  const double& operator()(size_t i) const { return storage_[i]; }

  size_t num_rows() const { return num_rows_; }

private:
  size_t              num_rows_;
  std::vector<double> storage_;
};

// From the original code
size_t FieldSolveMatrix(const Eigen::SparseMatrix<double> &A, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig, const double dx, const Vector &rho, Vector &phi, Vector &E_grid){
    size_t status = 0;
    // Initialize VectorXd's
    for (size_t ij = 0; ij < rhoEig.rows(); ij++){
        rhoEig[ij] = rho(ij); // A*phi = -rho
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
        phi(ij) = phiEig[ij];  
    }
    phi(phi.num_rows() - 1) = phi(0); // PBC

    // Calculate grid electric field 
    for (size_t ij = 0; ij < E_grid.num_rows(); ij++){
        if (ij == 0) // first-order central difference
            E_grid(ij) = -(phi(1) - phi(phi.num_rows()-2)) / (2.0 * dx);
        else if (ij == (E_grid.num_rows() - 1))
            E_grid(ij) = E_grid(0); // E_{Nx-1} = E_{0} is second PBC
            // E_grid(ij) = (phi(phi.num_rows() - 2) - phi(0)) / (2.0 * dx);
        else E_grid(ij) = -(phi(ij + 1) - phi(ij - 1)) / (2.0 * dx); 
    }
    return status;
}
#endif