#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cstddef>

// 11/21/23 - Needs to be refactored to work with PETSc
class Grid1d1v{
    public:
        Grid1d1v(size_t Nx) : 
        Nx_(Nx), x_grid(Nx, 0.0), rho_x(Nx, 0.0), E_x(Nx, 0.0)
        {}

        const size_t getNx() const { return Nx_; }

        const double RhoX(size_t j) const { return rho_x[j]; }
        const double EX(size_t j) const { return E_x[j]; }

        double& RhoX(size_t j) { return rho_x[j]; }
        double& EX(size_t j) { return E_x[j]; }

        const std::vector<double>& getRhoX() const { return rho_x; }
        const std::vector<double>& getEX() const { return E_x; }

        std::vector<double>& getRhoX() { return rho_x; }
        std::vector<double>& getEX() { return E_x; }

    private:
        size_t Nx_;
        std::vector<double> x_grid;
        std::vector<double> rho_x;
        std::vector<double> E_x; 
};
#endif