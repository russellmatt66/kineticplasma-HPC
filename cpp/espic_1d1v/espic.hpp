#ifndef ESPIC_HPP
#define ESPIC_HPP

#include<cmath>
#include<fstream>
#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>

#include "Vector.hpp"

using namespace std;

#define Q_particle -1.0

// Binary search - Ordered (sorted) grid
// Finds left gridpoint (tested)
size_t findParticle(const Vector &x_grid, const double particlePos){ 
    size_t low = 0;
    size_t high = x_grid.num_rows() - 1;
    size_t index;
    while (low <= high){
        index = floor((high + low) / 2); 
        if (x_grid(index) <= particlePos && x_grid(index + 1) >= particlePos) {
            return index;
        } else if (x_grid(index) > particlePos){
            high = index;
        } 
        else if (x_grid(index + 1) <= particlePos) { // better way than the below 'else if' is just to use '<=' here
            low = index + 1;
        } 
    }
    return -1; // not found
}

// trapezoidal integration for quasineutrality
double integrate(const Vector &data, const double dx){
    double sum = 0.0;
    for (size_t id = 0; id < data.num_rows()-1; id++){
        sum += 0.5 * (data(id) + data(id+1));
    }
    return dx * sum;
}

// Quasineutrality is fixed
// Current issue: What should particle_indices(ii) be for an out-of-bounds particle?
size_t ParticleWeight(const size_t W, Vector &particle_x, const Vector &x_grid, Vector &rho, Vector &particle_indices, double &Q_net)
{ 
    /* 
    @brief: Computes the charge density based on the distribution of the particles.     
    */
    size_t j_right, j_left, status = 0;
    double dist_left, dist_right;
    const double dx = (x_grid(x_grid.num_rows() - 1) - x_grid(0)) / (x_grid.num_rows() - 1); // a_0 = dx

    // no particles on grid = no charge density
    for (size_t ij = 0; ij < rho.num_rows(); ij++)
        rho(ij) = 0.0;

    // cout << "rho initialized" << endl;

    for (size_t ii = 0; ii < particle_x.num_rows(); ii++)
    {
        j_left = findParticle(x_grid, particle_x(ii)); // left gridpoint is found
        particle_indices(ii) = j_left;
        j_right = j_left + 1;
        dist_left = fabs(particle_x(ii) - x_grid(j_left));
        dist_right = fabs(x_grid(j_right) - particle_x(ii));  
        if (W == 0) 
        { // 0th-order weighting
            if (dist_right > dist_left) // closer to j_left
                rho(j_left) += Q_particle;
            else // must be closer to j_right, assigning "equals" case here is also fine, just a matter of where noise goes
                rho(j_right) += Q_particle;
        }
        else if (W == 1) 
        { // 1st-order weighting
            rho(j_left) += (Q_particle / dx) * dist_right;
            rho(j_right) += (Q_particle / dx) * dist_left;
        } 
    }
    
    rho(0) += rho(rho.num_rows()-1);
    rho(rho.num_rows()-1) = rho(0);

    // Quasineutrality
    // cout << "Establishing quasineutrality" << endl;
    Q_net = integrate(rho,dx); // this is negative charge
    for (size_t ij = 0; ij < rho.num_rows(); ij++) {
        rho(ij) += fabs(Q_net) / (x_grid(x_grid.num_rows() - 1) - x_grid(0));
    }
    Q_net = integrate(rho,dx);
    return status;
}

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

// 
size_t FieldSolveMatrix(const Eigen::SparseMatrix<double> &A, const double dx, const Vector &rho, Vector &phi, Vector &E_grid){
    size_t status = 0;
    Eigen::VectorXd rhoEig(A.rows()),phiEig(A.rows()); // Eigen interfaces for use with sparse solver
    // Initialize VectorXd's
    for (size_t ij = 0; ij < rhoEig.rows(); ij++){
        rhoEig[ij] = rho(ij); // A*phi = -rho
        phiEig[ij] = 0.0; // just initialize phi to 0 for simplicity
    }  
    // Boundary conditions - see writeup
    rhoEig[rhoEig.rows()-1] = 0.0;

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

// 
size_t ForceWeight(const size_t W, const Vector &particle_indices, const Vector &E_grid, const Vector &x_grid, const Vector &particle_x, Vector &particle_E){
    size_t j_right, j_left, status = 0;
    double dist_left, dist_right;
    const double dx = (x_grid(x_grid.num_rows() - 1) - x_grid(0)) / (x_grid.num_rows() - 1);

    for (size_t ii = 0; ii < particle_x.num_rows(); ii++){
        j_left = particle_indices(ii); 
        j_right = j_left + 1; // it's more efficient to just pass in the j-values from charge-weighting
        dist_left = fabs(particle_x(ii) - x_grid(j_left));
        dist_right = fabs(x_grid(j_right) - particle_x(ii));
        if (W == 0){ // 0th-order weighting
            if (dist_left < dist_right) particle_E(ii) = E_grid(j_left);
            else particle_E(ii) = E_grid(j_right);             
        }
        if (W == 1){ // 1st-order weighting
            particle_E(ii) = (dist_right / dx) * E_grid(j_left) + (dist_left / dx) * E_grid(j_right); 
            particle_E(ii) *= Q_particle;
        }
    }
    return status;
}

size_t ParticlePush(const Vector &particle_E, const double dt, const Vector &x_grid, Vector &particle_x, Vector &particle_v){
    size_t status = 0;
    double dx = (x_grid(x_grid.num_rows()-1)-x_grid(0))/x_grid.num_rows(), dist_right, dist_left; 
    for (size_t ii = 0; ii < particle_x.num_rows(); ii++){
        particle_v(ii) += dt * particle_E(ii); // compute v^{n+1/2}
        particle_x(ii) += dt * particle_v(ii); // use v^{n+1/2}
        // PBCs - update back into bounds
        if (particle_x(ii) > x_grid(x_grid.num_rows()-1)){ // oob to the right
            dist_right = fabs(particle_x(ii) - x_grid(x_grid.num_rows()-1)); // right of gridpoint
            particle_x(ii) = x_grid(0) + dist_right;
        }
        else if (particle_x(ii) < x_grid(0)) { // oob to the left
            dist_left = fabs(particle_x(ii) - x_grid(0)); // left of gridpoint
            particle_x(ii) = x_grid(x_grid.num_rows()-1) - dist_left;
        }
    }
    return status;
}

/* Check requirements specification for the way this needs to be handled */
size_t CollectParticleData(ofstream &phasespacefile, const char* filename, const double t, const size_t n, const Vector &particle_x, const Vector &particle_v, const Vector &particle_E, const Vector &particle_indices, const Vector &x_grid){
    size_t status = 0;
    // cout << "Collecting Particle Data" << endl;
    phasespacefile.open(filename, ios::app);
    for (size_t ii = 0; ii < particle_x.num_rows(); ii++){
        phasespacefile << t << "," << n << "," << ii << "," << particle_indices(ii) << "," << x_grid(particle_indices(ii)) << "," << particle_x(ii) << "," << particle_v(ii) << "," << particle_E(ii) << endl;
    }
    phasespacefile.close();
    // cout << "Particle Data Collected" << endl;
    return status;
}

size_t CollectGridData(ofstream &gridfile, const char* filename, const double t, const size_t n, const Vector &x_grid, const Vector &rho, const Vector &phi, const Vector &E_grid, const double Q_net){
    size_t status = 0;
    // cout << "Collecting Grid Data" << endl;
    gridfile.open(filename, ios::app); 
    for (size_t ij = 0; ij < x_grid.num_rows(); ij++){
        gridfile << t << "," << n << "," << ij << "," << x_grid(ij) << "," << rho(ij) << "," << phi(ij) << "," << E_grid(ij) << "," << Q_net << endl;
    }
    gridfile.close();
    // cout << "Grid Data Collected" << endl;
    return status;
}
// 
size_t CollectEnergyHistory(ofstream &energyfile, const char* filename, const double t, const size_t n, const double dt, const Vector &E_grid, const Vector &particle_v, const Vector &x_grid){
    size_t status = 0;
    // cout << "Collecting Energy History" << endl;
    energyfile.open(filename, ios::app); // append data
    double KEtotal = 0.0, Ftotal = 0.0, Etotal = 0.0; // system kinetic energy, (electric) field energy, total energy
    for (size_t ii = 0; ii < particle_v.num_rows(); ii++){ // sum KE
        KEtotal += 0.5 * pow(particle_v(ii),2);  
    }
    double dx = (x_grid(x_grid.num_rows()-1) - x_grid(0)) / (x_grid.num_rows()-1);
    for (size_t ij = 0; ij < E_grid.num_rows(); ij++){
        Ftotal += 0.5 * pow(E_grid(ij),2) * dx; // this is really large for some reason
    } 
    Etotal = KEtotal + Ftotal;
    energyfile << t << "," << n << "," << dt << "," << KEtotal << "," << Ftotal << "," << Etotal << endl;
    energyfile.close();
    // cout << "Energy History Collected" << endl;
    return status;
}

#endif