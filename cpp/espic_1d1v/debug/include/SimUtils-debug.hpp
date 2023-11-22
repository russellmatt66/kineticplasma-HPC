#ifndef UTILS_DEBUG_H
#define UTILS_DEBUG_H

#include <string>

#include "Grid-debug.hpp"
#include "Particles-debug.hpp"

using std::string;

// Initial Conditions
// Species
size_t ParticleInitialConditions(ParticleSpecies1d1v& PS, const double vprime, const string particleModeICs){
    size_t status = 0;

    return status;
}

// Grid
size_t GridInitialConditions(Grid1d1v& Grid, const string gridModeICs){

}

// Collect Data

#endif