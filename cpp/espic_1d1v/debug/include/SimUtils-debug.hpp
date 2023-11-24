#ifndef UTILS_DEBUG_H
#define UTILS_DEBUG_H
#define _USE_MATH_DEFINES

#include <string>
#include <cmath>

#include "Grid-debug.hpp"
#include "Particles-debug.hpp"

using std::string;

// Initial Conditions
size_t InitialConditions(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const double vprime, const string particleModeICs){
    size_t status = 0;
    double dx = Grid.getDX();
    status = ParticleICs(PS, Grid, vprime, dx, particleModeICs);
    return status;
}

// Species
// Grid is needed to get relevant attributes, e.g., L
size_t ParticleICs(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const double vprime, const double dx, const string particleModeICs){
    size_t status = 0;

    if (particleModeICs == "random"){
        status = ParticleICsRandom(PS, vprime);
    }
    else if (particleModeICs == "uniform"){
        status = ParticleICsUniform(PS, Grid, vprime, dx);
    }

    return status;
}

// Initialize a Maxwellian Population with a specific FWHM
size_t ParticleICsRandom(ParticleSpecies1d1v& PS, const double vprime){
    size_t status = 0;

    return status;
}

// Initialize a uniformly spaced population with a sinusoidal velocity perturbation
size_t ParticleICsUniform(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const double vprime){
    size_t status = 0;
    
    // Get all the necessary attributes 
    size_t N = PS.getParticleNum(), Nx = Grid.getNx();
    double L = Grid.getL();
    double k = 2.0 * M_PI / L, x_max = Grid.Xgrid(Nx-1), x_min = Grid.Xgrid(0); 
    double dx_p = (x_max - x_min) / N;
    double m = (1.0 / (N-1)) * (x_max - x_min - dx_p), b = x_min + dx_p / 2.0;

    // Initialize uniformly spaced population with sinusoidal velocity perturbation
    for (size_t ii = 0; ii < N; ii++){
        PS.ParticleX(ii) =  m * ii + b;
        PS.ParticleVx(ii) = vprime*sin(k*(PS.ParticleX(ii)));
    }

    return status;
}

// Grid
size_t GridICs(Grid1d1v& Grid, const string gridModeICs){
    size_t status = 0;
    /* Is this even necessary? */
    return status;
}

// Collect Data

#endif