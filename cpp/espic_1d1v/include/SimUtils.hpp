#ifndef UTILS_DEBUG_H
#define UTILS_DEBUG_H
#define _USE_MATH_DEFINES

#include <string>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <fstream>

#include "Grid.hpp"
#include "Particles.hpp"

using std::string;

// Initial Conditions
// Species
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

// ICs for two-stream instability
// Only one electron population b/c then there's no need to refactor further
size_t ParticleICsTwostream(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const double vprime){
    size_t status = 0;

    // Get all the necessary attributes
    size_t N = PS.getParticleNum(), Nx = Grid.getNx();
    double L = Grid.getL();
    double k = 2.0 * M_PI / L, x_max = Grid.Xgrid(Nx-1), x_min = Grid.Xgrid(0); 
    double dx_p = (x_max - x_min) / N;
    double m = (1.0 / (N-1)) * (x_max - x_min - dx_p), b = x_min + dx_p / 2.0;

    // Initialize counterstreaming electrons
    // Uniformly-spaced populations with a sinusoidal perturbation that is 90 deg out of phase b/w the two
    // Streaming right: i \in [0, N/2 - 1]
    // Streaming left: i \in [N/2, N - 1]
    size_t l_index; 
    k *= 2.0; // equivalent to changing wavelength of perturbation
    for (size_t ii = 0; ii < N / 2; ii++){
        // Streaming right
        PS.ParticleVx(ii) = vprime;
        PS.ParticleX(ii) =  m * ii + b;
        PS.ParticleX(ii) += 0.00001*sin(k * PS.ParticleX(ii));
        // Streaming left
        l_index = ii - 1 + N / 2;
        PS.ParticleVx(l_index) = -vprime;
        PS.ParticleX(l_index) = m * ii + b;
        PS.ParticleX(l_index) += 0.00001*sin(k * PS.ParticleX(l_index) + M_PI / 2.0);
    }
    return status;
}

// Initialize a Maxwellian Population with a specific FWHM
size_t ParticleICsRandom(ParticleSpecies1d1v& PS, const double vprime){
    size_t status = 0;
    /* TODO */
    return status;
}

// Wrapper for the specific particle configurations
// Grid is needed to get relevant attributes, e.g., L
size_t ParticleICs(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const double vprime, const double dx, const string particleModeICs){
    size_t status = 0;
    // size_t Nx = Grid.getNx();

    if (particleModeICs == "random"){
        status = ParticleICsRandom(PS, vprime);
    }
    else if (particleModeICs == "uniform"){
        status = ParticleICsUniform(PS, Grid, vprime);
    }
    else if (particleModeICs == "twostream"){
        status = ParticleICsTwostream(PS, Grid, vprime);
    }
    return status;
}

// Put it all togther
size_t InitialConditions(ParticleSpecies1d1v& PS, Grid1d1v& Grid, const double vprime, const string particleModeICs){
    size_t status = 0;
    double dx = Grid.getDX();
    status = ParticleICs(PS, Grid, vprime, dx, particleModeICs);
    return status;
}


// Grid
size_t GridICs(Grid1d1v& Grid, const string gridModeICs){
    size_t status = 0;
    /* Is this even necessary? */
    return status;
}



// Collect Data
size_t CollectData(const ParticleSpecies1d1v& PS, const Grid1d1v& Grid, const size_t timelevel){
    size_t status = 0;
    string fileString = "../data/Particles/Particles_" + std::to_string(timelevel) + ".csv";
    std::ofstream ParticleFile(fileString);

    size_t N = PS.getParticleNum();
    if (ParticleFile.is_open()){
        ParticleFile << "i,x_i,v_i,E_i" << std::endl;
        for (size_t i = 0; i < N; i++){
            /* Stream data to ParticleFile */
            ParticleFile << i << "," << PS.ParticleX(i) << "," << PS.ParticleVx(i) << "," << PS.ParticleEx(i) << std::endl;
        } 
    }

    fileString = "../data/Grid/Grid_" + std::to_string(timelevel) + ".csv";
    std::ofstream GridFile(fileString);
    size_t Nx = Grid.getNx();
    if (GridFile.is_open()){
        GridFile << "j,x_j,rho_j,phi_j,E_j" << std::endl;
        for (size_t j = 0; j < Nx; j++){
            /* Stream data to GridFile */
            GridFile << j << "," << Grid.Xgrid(j) << "," << Grid.RhoX(j) << "," << Grid.PhiX(j) << "," << Grid.EX(j) << std::endl;
        } 
    }

    return status;
}

size_t CollectEnergyHistory(std::ofstream& EnergyFile, const ParticleSpecies1d1v& Particles, const Grid1d1v& Grid, std::vector<double> E_sq, const size_t timelevel){
    size_t status = 0;
    double PE = 0, KE = 0, E = 0; // potential (electrostatic field) energy, particle kinetic energy, total energy

    // Compute field energy pointwise 
    size_t Nx = Grid.getNx();
    for (size_t j = 0; j < Nx; j++){
        E_sq[j] = 0.5 * Grid.EX(j) * Grid.EX(j);
    }

    // Trapezoidal integration of field energy
    double dx = Grid.getDX();
    for (size_t j = 0; j < Nx - 1; j++){
        PE += 0.5 * (E_sq[j] + E_sq[j+1]); 
    }
    PE *= dx;

    // Compute particle kinetic energy
    size_t N = Particles.getParticleNum();
    for (size_t i = 0; i < N; i++){
        KE += 0.5 * Particles.ParticleVx(i) * Particles.ParticleVx(i);
    }

    E = PE + KE;
    
    // Stream data to Energy file
    EnergyFile << timelevel << "," << KE << "," << PE << "," << E << std::endl;

    return status;
}

#endif