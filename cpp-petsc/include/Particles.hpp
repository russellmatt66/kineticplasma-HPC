#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <cstddef>

class ParticleSpecies1d1v{
    public:
        ParticleSpecies1d1v(size_t N, double Q) : 
        N_(N), Q_(Q), particle_x(N, 0.0), particle_vx(N, 0.0), particle_Ex(N, 0.0) 
        {}

        const size_t getParticleNum() const { return N_; }
        
        // Element-wise accessor methods
        const double ParticleX(size_t i) const { return particle_x[i]; }
        const double ParticleVx(size_t i) const { return particle_vx[i]; }
        const double ParticleEx(size_t i) const {return particle_Ex[i]; }

        double &ParticleX(size_t i) { return particle_x[i]; }
        double &ParticleVx(size_t i) { return particle_vx[i]; }
        double &ParticleEx(size_t i) { return particle_Ex[i]; }

        // Also need to be able to access whole data members
        const std::vector<double>& getParticleX() const { return particle_x; }
        const std::vector<double>& getParticleVx() const { return particle_vx; }
        const std::vector<double>& getParticleEx() const { return particle_Ex; }

        std::vector<double>& getParticleX() { return particle_x; }
        std::vector<double>& getParticleVx() { return particle_vx; }
        std::vector<double>& getParticleEx() { return particle_Ex; }

    private:
        size_t N_;
        double Q_; 
        std::vector<double> particle_x;
        std::vector<double> particle_vx;
        std::vector<double> particle_Ex; // x-dir acceleration
};
#endif