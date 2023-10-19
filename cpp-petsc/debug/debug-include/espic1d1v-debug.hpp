#ifndef ESPIC_HPP_DEBUG
#define ESPIC_HPP_DEBUG

#include<cmath>

#include "Grid-debug.hpp"
#include "Particles-debug.hpp"

// Binary search for findParticles() below
// Finds gridpoint to the left of where particle is (the cell)
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

// Computes ParticleSpecies1d1v.x_found as a side-effect - actually just a waste of time
// void findParticles(ParticleSpecies1d1v &PS, Grid1d1v &Grid){
//     size_t found_index;
//     size_t N = PS.getParticleNum();
//     for (size_t in = 0; in < N; in++){
//         found_index = findParticle(PS.ParticleX(in), Grid.getXgrid());
//         PS.XFound(in) = found_index;
//     }   
// }

// Weight the particles to the grid 
void ParticleWeight(size_t W, ParticleSpecies1d1v &PS, Grid1d1v &Grid){
    size_t j_left, j_right;
    double dist_left, dist_right;
    size_t Nx = Grid.getNx(), N = PS.getParticleNum();
    const double dx = (Grid.Xgrid(Nx - 1) - Grid.Xgrid(0)) / (Grid.Xgrid(Nx - 1));
    const double Q_particle = PS.getParticleQ();
    for (size_t ii = 0; ii < N;; ii+){
        j_left = findParticle(PS.ParticleX(ii), Grid.getXgrid());
        PS.XFound(ii) = j_left;
        j_right = j_left + 1;
        dist_left = fabs(PS.getParticleX(ii) - Grid.Xgrid(j_left));
        dist_right = fabs(Grid.Xgrid(j_right) - PS.getParticleX(ii));
        if (W == 0) // Zeroth-order weighting
        {
            if (dist_right > dist_left)
                Grid.RhoX(j_left) += Q_particle;
            else    
                Grid.RhoX(j_right) += Q_particle;
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

// Field solve


// Weight the force back to the particles

// Push the particles

// Collect data
#endif