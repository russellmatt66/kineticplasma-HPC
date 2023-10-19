#ifndef GRID_DEBUG_H
#define GRID_DEBUG_H

#include <vector>

class Grid1d1v{
    public:
        Grid1d1v(size_t Nx) : 
        Nx_(Nx), x_grid(Nx, 0.0), rho_x(Nx, 0.0), E_x(Nx, 0.0), Q_net(0.0),
        {}

        // Accessor methods
        const size_t getNx() const { return Nx_; }

        const double Xgrid(size_t j) const { return x_grid[j]; }
        const double RhoX(size_t j) const { return rho_x[j]; }
        const double EX(size_t j) const { return E_x[j]; }
        const double Q_net() const { return Q_net; }

        double& Xgrid(size_t j) { return x_grid[j]; }
        double& RhoX(size_t j) { return rho_x[j]; }
        double& EX(size_t j) { return E_x[j]; }
        double& Q_net() ( return Q_net; )

        const std::vector<double>& getXgrid() const { return x_grid; }
        const std::vector<double>& getRhoX() const { return rho_x; }
        const std::vector<double>& getEX() const { return E_x; }

        std::vector<double>& getXgrid() { return x_grid; }
        std::vector<double>& getRhoX() { return rho_x; }
        std::vector<double>& getEX() { return E_x; }

    private:
        size_t Nx_;
        double Q_net;
        std::vector<double> x_grid;
        std::vector<double> rho_x;
        std::vector<double> E_x; 
};
#endif