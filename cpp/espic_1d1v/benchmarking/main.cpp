/* 
Specifications found on Canvas
Idea here is to compute quantities with C++ and then plot them with Python. 
To accomplish this, the C++ code will write the necessary data to a set of CSV files.
This CSV files will then be opened by a Python program that extracts the data and plots it.
To fully automate this process, piping can be employed in a shell script. 

Project implements an electrostatic, 1D1V Particle-in-Cell (PIC) simulation, and studies the Leapfrog instability. 
The PIC method is a kind of kinetic plasma simulation that evolves the motion of charged superparticles who are subject to electromagnetic forces.
A self-consistent problem is solved whereby the state of the superparticles is used to determine the charge and currents on the grid.
These quantities are then used to obtain the electric and magnetic fields. In order to advance the state of the particles, the fields are 
weighted back to them in order to compute the accelerating force. Afterwards, the particles are pushed using a Leapfrog method which requires 
to be jumpstarted with a half-step backwards.  
*/ 
#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<random>
#include<chrono>

// Solved this issue with compiler flag -I path/to/eigen
#include <Eigen/Sparse>

#include "../espic.hpp"
#include "../Vector.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])
{
    size_t N = stoul(argv[1]); // number of particles
    size_t Nx = stoul(argv[2]); // number of gridpoints (# grid cells = Nx - 1)
    size_t Nt = stoul(argv[3]); // number of timesteps
    size_t W = stoul(argv[4]); // particle/force-weighting order (0th or 1st)
    double vprime = stod(argv[5]); // sets initial velocity for N = 2, velocity perturbation for N = 64

    auto start_total = high_resolution_clock::now();

    cout << "Number of particles " << N << endl;
    cout << "Number of gridpoints " << Nx << endl;
    cout << "Number of timesteps " << Nt << endl;
    cout << "Weighting order " << W << endl;

    // Create output files
    ofstream energy_data, grid_data, particle_data;
    
    // Initialize uniform 1D grid 
    Vector x_grid = Vector(Nx), rho = Vector(Nx), phi = Vector(Nx), E_grid = Vector(Nx);
    double x_min = -M_PI, x_max = M_PI;
    double dx = ( x_max - x_min ) / (x_grid.num_rows() - 1);
    cout << "Grid spacing is " << dx << endl;
    double Q_net;

    for (size_t ij = 0; ij < x_grid.num_rows(); ij++){
        x_grid(ij) = x_min + ij * dx; // uniform grid
        rho(ij) = 0.0;
        phi(ij) = 0.0;
        E_grid(ij) = 0.0;
    }

    // Initialize time
    cout << "Initializing time" << endl;
    double omega_p = sqrt(N/(x_max - x_min)); // sqrt(N/L)
    double tau_p = (2.0 * M_PI) / omega_p;
    double t = 0.0; // long doubles are a nightmare
    double dt = 0.01 * tau_p; // should this be tau_p or omega_p?
    cout << "Timestep is " << dt << endl;
    cout << "omega_p * dt is " << omega_p*dt << endl;
    size_t n = 0; // timelevel

    // Initialize population of (negative) charged particles in phase-space.
    Vector particle_x = Vector(N), particle_v = Vector(N), particle_E = Vector(N), particle_indices = Vector(N); 
    
    // Start random number generator
    random_device rd; 
    mt19937 gen(rd()); // Mersenne Twister engine
    
    uniform_real_distribution<double> uniformdist(x_min,x_max);

    normal_distribution<double> maxwellian(0.0, 1.0 / sqrt(2.0 * log(2.0))); // FWHM gives T/m

    if (N == 2){ // uniformly spaced        
        if (!vprime){ // v(0) = 0
            particle_x(0) = -M_PI / 4.0;
            particle_x(1) = M_PI / 4.0;
            particle_v(0) = 0.0;
            particle_v(1) = 0.0;
        } else if (vprime){ // v(0) = v'
            // double vprime = 0.01; // base this on a thermal velocity that gives FWHM of 2.0?
            particle_x(0) = -M_PI / 2.0;
            particle_x(1) = M_PI / 2.0;
            particle_v(0) = vprime;
            particle_v(1) = -vprime; 
        }
    } else if (N == 64){ // sinusoidal perturbation
        double k = 1.0;
        double dx_p = (x_max - x_min) / N;
        double m = (1.0 / (N-1)) * (x_max - x_min - dx_p), b = x_min + dx_p / 2.0;
        for (size_t ii = 0; ii < particle_x.num_rows(); ii++){
            // particle_x(ii) = x_min + ii * dx_particles;
            particle_x(ii) = m * ii + b;
            particle_v(ii) = vprime*sin(k*(particle_x(ii)));
        }
    } else {
        for (size_t ii = 0; ii < particle_x.num_rows(); ii++){
            particle_x(ii) = uniformdist(gen);
            particle_v(ii) =  maxwellian(gen);
        }
    }
    /*
    Try: 
    Handle the 'non-magic' particle numbers
    Default x: uniform
    Default v: maxwellian
    */
    cout << "Particles finished initializing" << endl;
    cout << "Preparing datafiles" << endl;
    char *energy_string, *grid_string, *particle_string;
    if (N == 2 && !vprime){
        energy_string = "pythoncode/energy_historyN2v0.csv";
        grid_string = "pythoncode/grid_dataN2v0.csv";
        particle_string = "pythoncode/particles_phasespaceN2v0.csv";        
    } else if (N == 2 && vprime){
        energy_string = "pythoncode/energy_historyN2v1.csv";
        grid_string = "pythoncode/grid_dataN2v1.csv";
        particle_string = "pythoncode/particles_phasespaceN2v1.csv";
    } else if (N == 64){
        energy_string = "pythoncode/energy_historyN64.csv";
        grid_string = "pythoncode/grid_dataN64.csv";
        particle_string = "pythoncode/particles_phasespaceN64.csv";
    } else {
        energy_string = "pythoncode/energy_historytest.csv";
        grid_string = "pythoncode/grid_datatest.csv";
        particle_string = "pythoncode/particles_phasespacetest.csv";
    }
    energy_data.open(energy_string);
    grid_data.open(grid_string);
    particle_data.open(particle_string);

    cout << "Preparing energy datafile" << endl;
    energy_data << "t,n,dt,Kinetic Energy,Electric Energy,Total Energy\n"; 
    energy_data.close();
    
    cout << "Preparing grid datafile" << endl; 
    grid_data << "t,n,j,x_grid,rho,phi,E_j,Q_net\n";
    grid_data.close();  

    cout << "Preparing phasespace datafile" << endl; 
    particle_data << "t,n,i,j_left,x_grid(j_left),position,velocity,E_i\n"; // i is particle number
    particle_data.close();


    // Weight the initial distribution of particles to the grid to obtain initial rho
    size_t routineFlag = ParticleWeight(W, particle_x, x_grid, rho, particle_indices, Q_net); // compute rho as a side-effect
    /*
    Process routineFlag
    */ 
    // Obtain initial potential and field on the grid
    // Construct sparse matrix expressing linear system
    Eigen::SparseMatrix<double> A(Nx-1,Nx-1);
    A.reserve(Eigen::VectorXi::Constant(Nx-1,3));
    routineFlag = BuildSparseLapl(A,dx);
    /*
    Process routineFlag 
    */
    // Given initial charge configuration, obtain initial potential
    routineFlag = FieldSolveMatrix(A,dx,rho,phi,E_grid);
    /*
    Process routineFlag 
    */
    // Weight the initial grid field to the particles
    routineFlag = ForceWeight(W,particle_indices,E_grid,x_grid,particle_x,particle_E);  
    /*
    Process routineFlag
    */
    // Half-step backwards in velocity, then push the particles
    for (size_t in = 0; in < particle_v.num_rows(); in++){ 
        particle_v(in) = particle_v(in) - 0.5*dt*particle_E(in); 
    }
    routineFlag = ParticlePush(particle_E,dt,x_grid,particle_x,particle_v);
    /*
    Process routineFlag
    */

    // PIC Algorithm
    auto start_pic = high_resolution_clock::now();
    for (size_t it = 1; it < Nt; it++){
        t += dt;
        // Charge Weight
        routineFlag = ParticleWeight(W, particle_x, x_grid, rho, particle_indices, Q_net);  
        // Field Solve
        routineFlag = FieldSolveMatrix(A,dx,rho,phi,E_grid);  
        // Force Weight
        routineFlag = ForceWeight(W,particle_indices,E_grid,x_grid,particle_x,particle_E);
        // Particle Push
        routineFlag = ParticlePush(particle_E,dt,x_grid,particle_x,particle_v);
    }
    auto stop_pic = high_resolution_clock::now();
    auto duration_pic = duration_cast<milliseconds>(stop_pic - start_pic);

    cout << Nt << " timesteps of PIC algorithm processed. " << Nx << " grid points, and " << N << " particles" << endl;
    cout << "Completed in " << duration_pic.count() << " milliseconds" << endl;
    
    // Program completion
    cout << "Program halting successfully. Closing datafiles."  << endl;
    if (energy_data.is_open())
        energy_data.close();
    if (grid_data.is_open())
        grid_data.close();
    if (particle_data.is_open())
        particle_data.close();
    cout << "Datafiles successfully closed, terminating execution." << endl;
    auto stop_total = high_resolution_clock::now();
    auto duration_total = duration_cast<milliseconds>(stop_total - start_total);
    cout << "Time taken by code: " << duration_total.count() << " milliseconds" << endl;
    cout << "wpdt is " << omega_p * dt << endl;
    // cout << A << endl;
    return 0;
}