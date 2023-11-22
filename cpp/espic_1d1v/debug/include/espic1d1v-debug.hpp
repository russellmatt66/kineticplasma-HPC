#ifndef ESPIC_HPP_DEBUG
#define ESPIC_HPP_DEBUG

#include<cmath>
#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>

#include "Grid-debug.hpp"
#include "Particles-debug.hpp"

/*
Binary search for ParticleWeight() below
Finds gridpoint to the left of where particle is, i.e., finds the cell the particle is in
*/ 
size_t findParticle(const double particlePos, const std::vector<double> &x_grid){
    size_t low = 0;
    size_t high = x_grid.size() - 1;
    size_t index; // gridpoint to the left of where particle is found, i.e., the cell
    while (low <= high){
        index = floor((high + low) / 2);
        if (x_grid[index] <= particlePos && x_grid[index + 1] >= particlePos) {
            return index;
        } else if (x_grid[index] > particlePos) {
            high = index;
        } else if (x_grid[index + 1] <= particlePos) {
            low = index + 1;
        }
    }
    return -1; // not found
}


/*
Weight the particles to the grid 
*/
void ParticleWeight(ParticleSpecies1d1v &PS, Grid1d1v &Grid, const size_t W, const size_t Nx, const size_t N, const double dx){
    size_t j_left, j_right;
    double dist_left, dist_right;
    const double Q_particle = PS.getParticleQ();
    for (size_t ii = 0; ii < N; ii++){
        j_left = findParticle(PS.ParticleX(ii), Grid.getXgrid());
        PS.XFound(ii) = j_left;
        j_right = j_left + 1;
        dist_left = fabs(PS.ParticleX(ii) - Grid.Xgrid(j_left));
        dist_right = fabs(Grid.Xgrid(j_right) - PS.ParticleX(ii));
        if (W == 0) // Zeroth-order weighting
        {
            if (dist_right > dist_left){
                Grid.RhoX(j_left) += Q_particle;
            }
            else{
                Grid.RhoX(j_right) += Q_particle;
            }
        }
        else if (W == 1) // First-order weighting
        {
            Grid.RhoX(j_left) += (Q_particle / dx) * dist_right;
            Grid.RhoX(j_right) += (Q_particle / dx) * dist_left;
        }    
    }

    // Periodic Boundary Conditions
    Grid.RhoX(0) += Grid.RhoX(Nx - 1);
    Grid.RhoX(Nx - 1) = Grid.RhoX(0);
} 

// Build Laplacian Finite Difference Stencil
size_t BuildSparseLapl(Eigen::SparseMatrix<double> &A, const double dx){
    size_t status = 0;
    size_t Ax = A.rows();
    
    /*
    Check if A.rows = A.cols
    */ 

    // Poissons Equation => Tridiagonal stencil
    for (size_t iA = 0; iA < Ax; iA++){
        if (iA > 0) A.insert(iA,iA-1) = 1.0;
        A.insert(iA,iA) = -2.0;
        if (iA < Ax - 1) A.insert(iA,iA+1) = 1.0; 
    }

    // Gauge
    for (size_t iA = 0; iA < Ax; iA++){
        A.coeffRef(Ax-1,iA) = 1.0;
    }

    // Periodic Boundaries
    A.coeffRef(0,Ax-1) = 1.0;
    
    A = (-1.0 / pow(dx,2)) * A;

    A.finalize();

    return status;
}

// Field Solve
size_t FieldSolveMatrix(const Eigen::SparseMatrix<double>& A, Grid1d1v& Grid, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig, const double dx, const size_t Nx){
    size_t status = 0;

    // Initialize VectorXd's
    for (size_t ij = 0; ij < Nx; ij++){
        rhoEig[ij] = Grid.RhoX(ij); 
        phiEig[ij] = 0.0; // just initialize phi to 0 for simplicity
    }  
    
    // Boundary conditions - see writeup
    rhoEig[Nx-1] = 0.0;

    // Solve A*phi = -rho for phi
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    phiEig = solver.solve(rhoEig);
    double residual = (A * phiEig - rhoEig).norm();
    // cout << "LU residual is " << residual << endl;

    // Copy out
    for (size_t ij = 0; ij < Nx; ij++){
        Grid.PhiX(ij) = phiEig[ij];  
    }

    // PBC
    Grid.PhiX(Nx - 1) = Grid.PhiX(0); 

    // Calculate grid electric field - first-order central difference
    for (size_t ij = 0; ij < Nx; ij++){
        if (ij == 0) {
            Grid.EX(ij) = -(Grid.PhiX(1) - Grid.PhiX(Nx - 2)) / (2.0 * dx);
        }
        else if (ij == (Nx - 1)){
            Grid.EX(ij) = Grid.EX(0); // E_{Nx-1} = E_{0} is second PBC
            // E_grid(ij) = (phi(phi.num_rows() - 2) - phi(0)) / (2.0 * dx);
        }
        else {
            Grid.EX(ij) = -(Grid.PhiX(ij + 1) - Grid.PhiX(ij - 1)) / (2.0 * dx); 
        }
    }
    return status;
}

// Weight the force back to the particles
size_t ForceWeight(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const size_t N, const size_t W, const double dx){
    size_t j_right, j_left, status = 0;
    double dist_left, dist_right;

    for (size_t ii = 0; ii < N; ii++){
        j_left = PS.XFound(ii); 
        j_right = j_left + 1; 
        dist_left = fabs(PS.ParticleX(ii) - Grid.Xgrid(j_left));
        dist_right = fabs(Grid.Xgrid(j_right) - PS.ParticleX(ii));
        if (W == 0){ // 0th-order weighting
            if (dist_left < dist_right) {
                PS.ParticleEx(ii) = Grid.EX(j_left);
            } 
            else {
                PS.ParticleEx(ii) = Grid.EX(j_right);             
            }
        }
        if (W == 1){ // 1st-order weighting
            PS.ParticleEx(ii) = (dist_right / dx) * Grid.EX(j_left) + (dist_left / dx) * Grid.EX(j_right); 
            PS.ParticleEx(ii) *= PS.getParticleQ();
        }
    }
    return status;
}

// Push the particles
// size_t ParticlePush(const Vector &particle_E, const double dt, const Vector &x_grid, Vector &particle_x, Vector &particle_v)
size_t ParticlePush(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const size_t N, const size_t Nx, const double dt, const double dx){
    size_t status = 0;
    double dist_right, dist_left;
    for (size_t ii = 0; ii < N; ii++){
        PS.ParticleVx(ii) += dt * PS.ParticleEx(ii); // compute v^{n+1/2}
        PS.ParticleX(ii) += dt * PS.ParticleVx(ii); // use v^{n+1/2} to compute x^{n+1}
        // PBCs - update back into bounds
        if (PS.ParticleX(ii) > Grid.Xgrid(Nx-1)){ // oob to the right
            dist_right = fabs(PS.ParticleX(ii) - Grid.Xgrid(Nx-1)); // distance overshot the right of gridpoint
            PS.ParticleX(ii) = Grid.Xgrid(0) + dist_right;
        }
        else if (PS.ParticleX(ii) < Grid.Xgrid(0)) { // oob to the left
            dist_left = fabs(PS.ParticleX(ii) - Grid.Xgrid(0)); // distance overshot the left of gridpoint
            PS.ParticleX(ii) = Grid.Xgrid(Nx-1) - dist_left;
        }
    }
    return status;
}
#endif