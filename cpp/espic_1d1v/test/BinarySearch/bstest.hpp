#ifndef BSTEST_HPP
#define BSTEST_HPP

#include "../../include/Grid.hpp"
#include "../../include/Particles.hpp"
#include "../../include/espic1d1v.hpp"

#include <fstream>
#include <vector>
#include <random>

// Validate results
size_t LinearSearch( const double ParticlePos, Grid1d1v& Grid){
    for (size_t j = 0; j < Grid.getNx()-1; j++){
        if (ParticlePos >= Grid.Xgrid(j) && ParticlePos < Grid.Xgrid(j+1)){
            return j;
        }
    }
    return -1; // not found 
}

// Case 1: Single particle in each cell
std::vector<bool> Cell(std::ofstream& ResultsLog, Grid1d1v& Grid, ParticleSpecies1d1v& Particles, const double x_min, const double x_max){
    size_t N = Particles.getParticleNum();
    double dx = Grid.getDX();
    double m = (x_max - x_min - dx) / (N-1), b = x_min + dx / 2.0;
    
    // Initialize single particle in each cell
    for (size_t i = 0; i < N; i++){
        Particles.ParticleX(i) = m * i + b;
    }

    size_t lsidx, bsidx;
    std::vector<bool> foundSame(N);
    for (size_t i = 0; i < N; i++){
        ResultsLog << "Searching for Particle " << i << " at x = " << Particles.ParticleX(i) << std::endl;
        lsidx = LinearSearch(Particles.ParticleX(i), Grid);
        bsidx = findParticle(Particles.ParticleX(i), Grid.getXgrid()); // binary search
        ResultsLog << "Linear Search found Particle " << i << " in Cell " << lsidx << std::endl;
        ResultsLog << "Binary Search found Particle " << i << " in Cell " << bsidx << std::endl;
        foundSame[i] = (lsidx == bsidx);
    }
    ResultsLog << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    return foundSame;
}

// Case 2: Randomly distributed particles
std::vector<bool> randomParticles(std::ofstream& ResultsLog, Grid1d1v& Grid, ParticleSpecies1d1v& Particles, const std::uniform_real_distribution<double> ufmrl_dist, const std::mt19937 gen){
    size_t N = Particles.getParticleNum();
    for (size_t i = 0; i < N; i++){
        Particles.ParticleX(i) = ufmrl_dist(gen);
    }

    size_t lsidx, bsidx; 
    std::vector<bool> foundSame(N);
    for (size_t i = 0; i < N; i++){
        ResultsLog << "Searching for Particle " << i << " at x = " << Particles.ParticleX(i) << std::endl;
        lsidx = LinearSearch(Particles.ParticleX(i), Grid);
        bsidx = findParticle(Particles.ParticleX(i), Grid.getXgrid()); // binary search
        ResultsLog << "Linear Search found Particle " << i << " in Cell " << lsidx << std::endl;
        ResultsLog << "Binary Search found Particle " << i << " in Cell " << bsidx << std::endl;
        foundSame[i] = (lsidx == bsidx);
    }
    ResultsLog << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    return foundSame;
}

// Case 3: Single particle on each left gridpoint
std::vector<bool> leftGridpoint(std::ofstream& ResultsLog, Grid1d1v& Grid, ParticleSpecies1d1v& Particles){
    size_t N = Particles.getParticleNum();

    for (size_t i = 0; i < N-1; i++){
        Particles.ParticleX(i) = Grid.Xgrid(i);
    }

    size_t lsidx, bsidx; 
    std::vector<bool> foundSame(N);
    for (size_t i = 0; i < N; i++){
        ResultsLog << "Searching for Particle " << i << " at x = " << Particles.ParticleX(i) << std::endl;
        lsidx = LinearSearch(Particles.ParticleX(i), Grid);
        bsidx = findParticle(Particles.ParticleX(i), Grid.getXgrid()); // binary search
        ResultsLog << "Linear Search found Particle " << i << " in Cell " << lsidx << std::endl;
        ResultsLog << "Binary Search found Particle " << i << " in Cell " << bsidx << std::endl;
        foundSame[i] = (lsidx == bsidx);
    }
    ResultsLog << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    return foundSame;
}

// Case 4: Single particle on each right gridpoint
std::vector<bool> rightGridpoint(std::ofstream& ResultsLog, Grid1d1v& Grid, ParticleSpecies1d1v& Particles){
    size_t N = Particles.getParticleNum();

    for (size_t i = 0; i < N-1; i++){
        Particles.ParticleX(i) = Grid.Xgrid(i+1);
    }

    size_t lsidx, bsidx; 
    std::vector<bool> foundSame(N);
    for (size_t i = 0; i < N; i++){
        ResultsLog << "Searching for Particle " << i << " at x = " << Particles.ParticleX(i) << std::endl;
        lsidx = LinearSearch(Particles.ParticleX(i), Grid);
        bsidx = findParticle(Particles.ParticleX(i), Grid.getXgrid()); // binary search
        ResultsLog << "Linear Search found Particle " << i << " in Cell " << lsidx << std::endl;
        ResultsLog << "Binary Search found Particle " << i << " in Cell " << bsidx << std::endl;
        foundSame[i] = (lsidx == bsidx);
    }
    ResultsLog << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    return foundSame;
}
#endif