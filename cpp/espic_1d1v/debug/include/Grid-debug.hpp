#ifndef GRID_DEBUG_H
#define GRID_DEBUG_H

#include <vector>

class Grid1d1v{
    public:
        Grid1d1v(size_t Nx, double x_min, double x_max) : 
        Nx_(Nx), Qnet_(0.0), x_grid(Nx, 0.0), rho_x(Nx, 0.0), phi_x(Nx, 0.0), E_x(Nx, 0.0)
        {
            // uniformly spaced grid
            double dx = (x_max - x_min) / (Nx - 1);
            for (size_t ij = 0; ij < Nx_; ij++){
                x_grid[ij] = x_min + ij * dx; 
            }
        }

        // Accessor methods
        const size_t getNx() const { return Nx_; }
        const size_t getQnet() const { return Qnet_; }
        const double getL() const { return x_grid[Nx_ - 1] - x_grid[0]; }
        const double getDX() const { return getL() / (Nx_ - 1); }

        const double Xgrid(size_t j) const { return x_grid[j]; }
        const double RhoX(size_t j) const { return rho_x[j]; }
        const double PhiX(size_t j) const { return phi_x[j]; }
        const double EX(size_t j) const { return E_x[j]; }

        double& Xgrid(size_t j) { return x_grid[j]; }
        double& RhoX(size_t j) { return rho_x[j]; }
        double& PhiX(size_t j) { return phi_x[j]; }
        double& EX(size_t j) { return E_x[j]; }

        const std::vector<double>& getXgrid() const { return x_grid; }
        const std::vector<double>& getRhoX() const { return rho_x; }
        const std::vector<double>& getEX() const { return E_x; }

        std::vector<double>& getXgrid() { return x_grid; }
        std::vector<double>& getRhoX() { return rho_x; }
        std::vector<double>& getEX() { return E_x; }

        // Zero out the charge density
        // Needed for ParticleWeighting() to not accumulate excess charge
        void ZeroOutRho() {
            for (size_t ij = 0; ij < Nx_; ij++){
                rho_x[ij] = 0.0;
            }
        }

        // Calculate net charge on grid
        // Trapezoidal integratation of rho_x across x_grid
        double calculate_Qnet() { 
            double dx = getDX(), sum = 0.0;
            for (size_t ij = 0; ij < Nx_-1; ij++){
                sum += 0.5 * (rho_x[ij] + rho_x[ij+1]);
            }
            Qnet_ = dx * sum;
            return Qnet_; 
        }

        // Neutralize grid with uniform, positive background
        void UniformPositiveBackground() {
            double L = getL(), totalQ = fabs(Qnet_);
            for (size_t ij = 0; ij < Nx_; ij++){
                rho_x[ij] += totalQ / L; 
            }
        }

    private:
        size_t Nx_;
        double Qnet_;
        std::vector<double> x_grid;
        std::vector<double> rho_x;
        std::vector<double> phi_x;
        std::vector<double> E_x; 
};
#endif