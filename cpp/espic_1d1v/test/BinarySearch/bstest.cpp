#define _USE_MATH_DEFINES

#include "../../include/Grid.hpp"
#include "../../include/Particles.hpp"
#include "bstest.hpp"

#include <fstream>
#include <random>

// Driver code for testing Binary Search
int main(int argc, char *argv[]){
    std::ofstream testlog("testsPassed.log"), resultslog("data.log");

    size_t Nx = std::stoi(argv[2]); // gonna test from 32 to 8192
    double x_min = -M_PI, x_max = M_PI;
    Grid1d1v testGrid(Nx,x_min,x_max);
    
    size_t N = std::stoi(argv[1]); // gonna test 64, and 1024. Testing involves linear search so keep small. 
    ParticleSpecies1d1v testParticles(N,-1.0);

    // Run different cases and process output
    // Cell
    testlog << "Beginning Cell test." << std::endl; 
    std::vector<bool> foundSame = Cell(resultslog, testGrid, testParticles, x_min, x_max);
    bool passed = true;
    for (int ib = 0; ib < foundSame.size(); ib++){
        if (!foundSame[ib]){
            passed = false;
        }
    }
    testlog << "Cell: " << passed << std::endl;

    // Random
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> ufmrl_dist(x_min, x_max);

    testlog << "Beginning Random test." << std::endl; 
    foundSame = randomParticles(resultslog, testGrid, testParticles, ufmrl_dist, gen);
    bool passed = true;
    for (int ib = 0; ib < foundSame.size(); ib++){
        if (!foundSame[ib]){
            passed = false;
        }
    }
    testlog << "Random: " << passed << std::endl;

    // Left gridpoint
    testlog << "Beginning leftGridpoint test." << std::endl; 
    ParticleSpecies1d1v testParticles(Nx - 1,-1.0); // need N = Nx - 1
    foundSame = leftGridpoint(resultslog, testGrid, testParticles);
    bool passed = true;
    for (int ib = 0; ib < foundSame.size(); ib++){
        if (!foundSame[ib]){
            passed = false;
        }
    }
    testlog << "leftGridpoint: " << passed << std::endl;

    // Right gridpoint - N = Nx
    testlog << "Beginning rightGridpoint test." << std::endl; 
    foundSame = rightGridpoint(resultslog, testGrid, testParticles);
    bool passed = true;
    for (int ib = 0; ib < foundSame.size(); ib++){
        if (!foundSame[ib]){
            passed = false;
        }
    }
    testlog << "rightGridpoint: " << passed << std::endl;

    return 0;
}