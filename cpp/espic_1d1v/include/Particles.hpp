#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>

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
#endif