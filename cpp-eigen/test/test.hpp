/* 
This is the test library.
It will define functions that:
(1) Aid in the testing of the project components 
*/ 
#ifndef TEST_HPP
#define TEST_HPP

#include<iostream>
#include<fstream>
#include<cmath>

#include "C:/ProgramData/eigen-3.4.0/Eigen/Sparse"
#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Eigenvalues>

#include "../Vector.hpp"

using namespace std;

#define Q_particle -1.0
#define eps_0 1.0
#define q_over_m -1.0
// const double Q_particle = -1.0; // unit charge
// const double eps_0 = 1.0;
// const double q_over_m = -1.0;

// Binary search - bugs should be fixed, truth conditional and logic 
size_t findParticleTest(const Vector &x_grid, const double particlePos){ // bug should be fixed
    cout << "inside Binary Search test function" << endl;
    size_t low = 0;
    size_t high = x_grid.num_rows() - 1;
    size_t index;
    size_t count = 0;
    while (low <= high && count < 10){
        count+=1; 
        index = floor((high + low) / 2); 
        cout << "Low is " << low << endl;
        cout << "High is " << high << endl;
        cout << "Current guess index is " << index << endl; 
        cout << "Particle location is " << particlePos << endl;
        cout << "x_grid(index) is " << x_grid(index) << endl;
        cout << "x_grid(index+1) is " << x_grid(index+1) << endl;
        if (x_grid(index) <= particlePos && x_grid(index + 1) > particlePos) {
            return index;
        } else if (x_grid(index) > particlePos){
            high = index;
        } 
        else if (x_grid(index + 1) < particlePos) {
            low = index + 1;
        } 
    }
    return -1; // not found
    // size_t low = 0;
    // size_t high = x_grid.num_rows() - 1;
    // size_t index;
    // while (low <= high){
    //     index = floor((high + low) / 2);
    //     cout << "Low is " << low << endl;
    //     cout << "High is " << high << endl;
    //     cout << "Current guess index is " << index << endl; 
    //     cout << "Particle location is " << particlePos << endl;
    //     cout << "x_grid(index) is " << x_grid(index) << endl;
    //     cout << "x_grid(index+1) is " << x_grid(index+1) << endl;
    //     if (x_grid(index) <= particlePos && x_grid(index + 1) > particlePos) { // bug was here - inside truth conditional
    //         return index;
    //     } else if (x_grid(index) > particlePos){
    //         high = index;
    //     } 
    //     else if (x_grid(index + 1) < particlePos) {
    //         low = index + 1;
    //     } 
    // }
    // return -1; // not found
}

double integrateTest(const Vector &data, const double dx){
    double sum = 0.0;
    for (size_t id = 0; id < data.num_rows()-1; id++){
        sum += 0.5 * (data(id) + data(id+1));
    }
    return dx * sum;
}

size_t findParticleRaster(const Vector &x_grid, const double particlePos){ // slow 
    for (size_t ix = 0; ix < x_grid.num_rows(); ix++){
        if (particlePos >= x_grid(ix) && particlePos < x_grid(ix+1)) return ix; // left gridpoint
    }
    return -1; // not found
}

size_t ParticleWeightTest(const size_t W, const Vector &particle_x, const Vector &x_grid, Vector &rho, double &Q_net)
{ 
    /* 
    @brief: Computes the charge density based on the distribution of the particles.     
    */
    cout << "Inside ParticleWeightTest " << endl;
    size_t j_right, j_left, status = 0;
    const double dx = (x_grid(x_grid.num_rows() - 1) - x_grid(0)) / (x_grid.num_rows() - 1); // a_0 = dx

    // if rho(ij) != 0 for all ij
    // then rho(ij) = 0 for all ij
    for (size_t ij = 0; ij < rho.num_rows(); ij++)
        rho(ij) = 0.0;

    for (size_t ii = 0; ii < particle_x.num_rows(); ii++)
    {
        double particlePos = particle_x(ii);
        j_left = findParticleTest(x_grid, particlePos); // left gridpoint is found
        j_right = j_left + 1;
        if (W == 0) { // 0th-order weighting
            if (fabs(particlePos - x_grid(j_left)) > fabs(particlePos - x_grid(j_right))) // closer to j_right
                rho(j_right) = Q_particle / dx;
            else // must be closer to j_left, assigning "equals" case to j_left is also fine, just a matter of where noise goes
                rho(j_left) = Q_particle / dx;
        }
        else if (W == 1) { // 1st-order weighting
            rho(j_right) = (Q_particle / dx) * fabs(particlePos - x_grid(j_left));
            rho(j_left) = (Q_particle / dx) * fabs(x_grid(j_right) - particlePos);
        } 
        
    }
    // Quasineutrality
    cout << "Establishing quasineutrality" << endl;
    Q_net = integrateTest(rho,dx); // this is negative charge
    for (size_t ij = 0; ij < rho.num_rows(); ij++) {
        // rho(ij) += fabs(Q_net) * (particle_x.num_rows() / (x_grid(x_grid.num_rows() - 1) - x_grid(0)));
        rho(ij) += fabs(Q_net) / (x_grid(x_grid.num_rows() - 1) - x_grid(0));
    }
    Q_net = integrateTest(rho,dx);
   return status;
}

// check that no charge = no field
size_t FieldSolveTest(){
    size_t status = 1;

    return status;
}
#endif