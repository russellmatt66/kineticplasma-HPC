#ifndef UTILS_DEBUG_H
#define UTILS_DEBUG_H

#include <string>

#include "Grid-debug.hpp"
#include "Particles-debug.hpp"

using std::string;

// Initial Conditions
size_t InitialConditions(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const double vprime, const string particleModeICs){
    size_t status = 0;

    return status;
}

// Species
size_t ParticleICs(ParticleSpecies1d1v& PS, const double vprime, const string particleModeICs){
    size_t status = 0;

    if (particleModeICs == "random"){
        status = ParticleICsRandom();
    }
    else if (particleModeICs == "uniform"){
        status = ParticleICsUniform();
    }

    return status;
}

// Initialize a Maxwellian Population with a specific FWHM
size_t ParticleICsRandom(){
    size_t status = 0;

    return status;
}

// Initialize a uniformly spaced population
size_t ParticleICsUniform(){
    size_t status = 0;

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