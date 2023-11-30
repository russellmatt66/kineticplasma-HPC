#define _USE_MATH_DEFINES

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <variant>
#include <string>
#include <cstddef>

#include <Eigen/Sparse>

#include "../include/Particles.hpp"
#include "../include/Grid.hpp"
#include "../include/espic1d1v.hpp"
#include "../include/SimUtils.hpp"

using std::string;
using ParameterValue = std::variant<size_t, double, string>;

std::unordered_map<string, ParameterValue> parseInputFile(const string filename);

int main(){
    size_t routineFlag;
    std::ofstream simlog;
    simlog.open("espic1d1v.log");
    simlog << "Beginning simulation" << std::endl;
    
    // Parse input file
    simlog << "Parsing input file" << std::endl;
    std::unordered_map<string, ParameterValue> inputParameters = parseInputFile("espic1d1v.inp");
    size_t N = std::get<size_t>(inputParameters["N"]);
    size_t Nx = std::get<size_t>(inputParameters["Nx"]);
    size_t Nt = std::get<size_t>(inputParameters["Nt"]);
    size_t W = std::get<size_t>(inputParameters["W"]);
    double vprime = std::get<double>(inputParameters["vprime"]);

    // Initialize simulation: particle species, and uniform grid
    simlog << "Constructing particles" << std::endl;
    ParticleSpecies1d1v electrons(N,-1.0);
    simlog << "Constructing grid" << std::endl;
    Grid1d1v Grid(Nx, -M_PI, M_PI);  
    double dx = Grid.getDX();
    double L = Grid.getL();


    string particleICs = std::get<string>(inputParameters["particleICs"]);
    simlog << "Writing Initial Conditions" << std::endl;
    routineFlag = InitialConditions(electrons,Grid,vprime,particleICs);

    // Initialize time
    simlog << "Initializing time variables" << std::endl;

    double omega_p = sqrt(N/L);
    double t = 0.0, tau_p = (2.0 * M_PI) / omega_p;
    double dt = std::get<double>(inputParameters["dtcoeff"]) * tau_p;
    
    /*
    Initial PIC step
    */ 
    simlog << "Performing initial PIC step" << std::endl;  
    routineFlag = ParticleWeight(electrons,Grid,W,Nx,N,dx);

    // Build sparse matrix for electrostatic field solve
    // A is (Nx-1)x(Nx-1) because of a persistent bug when it was (Nx)x(Nx)
    Eigen::SparseMatrix<double> A(Nx-1,Nx-1); // Periodic Boundary Conditions so don't need last line 
    A.reserve(Eigen::VectorXi::Constant(Nx-1,3)); // Poisson's equation so triangular
    routineFlag = BuildSparseLapl(A,dx);

    // Electrostatic field solve with a sparse matrix
    Eigen::VectorXd rhoEig(A.rows()), phiEig(A.rows()); // Eigen needs its own containers for sparse solver
    routineFlag = FieldSolveMatrix(A,Grid,rhoEig,phiEig,dx,Nx);

    // Weight grid electric field to particles, and then push 
    routineFlag = ForceWeight(electrons,Grid,N,W,dx); 
    routineFlag = CollectData(electrons,Grid,0); // Collect the data before pushing
    for (size_t ii = 0; ii < N; ii++){ // half-step backwards to start Leapfrog scheme
        electrons.ParticleVx(ii) = electrons.ParticleVx(ii) - 0.5*dt*electrons.ParticleEx(ii);
    }
    routineFlag = ParticlePush(electrons,Grid,N,Nx,dt,dx);



    // PIC Loop
    for (size_t it = 1; it < Nt; it++){
        simlog << "Starting timestep " << it << std::endl;
        t += dt;
        Grid.ZeroOutRho(); // Don't want to accumulate excess charge density between timesteps
        routineFlag = ParticleWeight(electrons,Grid,W,Nx,N,dx);
        routineFlag = FieldSolveMatrix(A,Grid,rhoEig,phiEig,dx,Nx);
        routineFlag = ForceWeight(electrons,Grid,N,W,dx); 
        routineFlag = CollectData(electrons,Grid,it);
        routineFlag = ParticlePush(electrons,Grid,N,Nx,dt,dx);
        simlog << "Timestep " << it << " complete" << std::endl;
    }

    simlog << "Closing log" << std::endl;
    simlog.close();
    return 0;
}

std::unordered_map<string, ParameterValue> parseInputFile(const string filename){
    std::unordered_map<string, ParameterValue> parameters;
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file: " << filename << std::endl;
        return parameters;
    }

    string line, paramName, paramValueStr;

    while(getline(inputFile, line)){
        size_t delimiterPos = line.find("=");
        if (delimiterPos != string::npos){
            string paramName = line.substr(0,delimiterPos);
            string paramValueStr = line.substr(delimiterPos + 1);
            if (paramName == "N" || paramName == "Nt" || paramName == "Nx" || paramName == "W"){
                size_t paramValue = std::stoul(paramValueStr);
                parameters[paramName] = paramValue;
            } 
            else if (paramName == "vprime" || paramName == "dtcoeff"){
                double paramValue = std::stod(paramValueStr);
                parameters[paramName] = paramValue;
            }
            else if (paramName == "particleICs"){
                parameters[paramName] = paramValueStr;
            }
        }
    }
    return parameters;
}