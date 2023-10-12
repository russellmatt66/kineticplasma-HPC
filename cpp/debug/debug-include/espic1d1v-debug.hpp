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

// Computes ParticleSpecies1d1v.x_found as a side-effect
void findParticles(ParticleSpecies1d1v &PS, Grid1d1v &Grid){
    size_t temp;
    size_t N = PS.getParticleNum();
    for (size_t in = 0; in < N; in++){
        temp = findParticle(PS.ParticleX(in), Grid.getXgrid());
        PS.XFound(in) = temp;
    }
}

#endif