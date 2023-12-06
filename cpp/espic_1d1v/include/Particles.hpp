#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <cstddef>
#include <string>

class ParticleSpecies1d1v{
    public:
        ParticleSpecies1d1v(size_t N, double Q) : 
        N_(N), Q_(Q), particle_x(N, 0.0), particle_vx(N, 0.0), particle_Ex(N, 0.0),
        x_found(N, 0) 
        {}

        const size_t getParticleNum() const { return N_; }
        const double getParticleQ() const { return Q_; }
        
        // Element-wise accessor methods
        const double ParticleX(size_t i) const { return particle_x[i]; }
        const double ParticleVx(size_t i) const { return particle_vx[i]; }
        const double ParticleEx(size_t i) const {return particle_Ex[i]; }
        const size_t XFound(size_t i) const { return x_found[i]; }

        double &ParticleX(size_t i) { return particle_x[i]; }
        double &ParticleVx(size_t i) { return particle_vx[i]; }
        double &ParticleEx(size_t i) { return particle_Ex[i]; }
        size_t &XFound(size_t i) { return x_found[i]; }

        // Also need to be able to access whole data members
        const std::vector<double>& getParticleX() const { return particle_x; }
        const std::vector<double>& getParticleVx() const { return particle_vx; }
        const std::vector<double>& getParticleEx() const { return particle_Ex; }
        const std::vector<size_t>& getXFound() const { return x_found; }

        std::vector<double>& getParticleX() { return particle_x; }
        std::vector<double>& getParticleVx() { return particle_vx; }
        std::vector<double>& getParticleEx() { return particle_Ex; }
        std::vector<size_t>& getXFound() { return x_found; }

    private:
        size_t N_;
        double Q_; 
        std::vector<double> particle_x;
        std::vector<double> particle_vx;
        std::vector<double> particle_Ex; // x-dir acceleration
        std::vector<size_t> x_found; // location of particles
};

// Class to collect all the particles together
class ParticleList{
    public:
        ParticleList(size_t N) : 
        N_(N), particles_x(N, 0.0), particles_vx(N, 0.0), particles_Ex(N, 0.0), x_found(N, 0), species_id(N, "") 
        {}

        // Boilerplate
        const size_t getParticleNum() const { return N_; }

        // Element-wise accessor methods
        const double ParticleX(size_t i) const { return particles_x[i]; }
        const double ParticleVx(size_t i) const { return particles_vx[i]; }
        const double ParticleEx(size_t i) const {return particles_Ex[i]; }
        const size_t XFound(size_t i) const { return x_found[i]; }

        double &ParticleX(size_t i) { return particles_x[i]; }
        double &ParticleVx(size_t i) { return particles_vx[i]; }
        double &ParticleEx(size_t i) { return particles_Ex[i]; }
        size_t &XFound(size_t i) { return x_found[i]; }

        // Also need to be able to access whole data members
        const std::vector<double>& getParticleX() const { return particles_x; }
        const std::vector<double>& getParticleVx() const { return particles_vx; }
        const std::vector<double>& getParticleEx() const { return particles_Ex; }
        const std::vector<size_t>& getXFound() const { return x_found; }

        std::vector<double>& getParticleX() { return particles_x; }
        std::vector<double>& getParticleVx() { return particles_vx; }
        std::vector<double>& getParticleEx() { return particles_Ex; }
        std::vector<size_t>& getXFound() { return x_found; }

    private:
        size_t N_;
        std::vector<double> particles_x;
        std::vector<double> particles_vx;
        std::vector<double> particles_Ex;
        std::vector<size_t> x_found;
        std::vector<std::string> species_id;
};
#endif