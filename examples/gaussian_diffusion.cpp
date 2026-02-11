#include "heat_solver.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>
#include <vector>

/**
 * @brief Example: Gaussian pulse diffusion
 * 
 * This example demonstrates the diffusion of a Gaussian pulse
 * using the 1D heat equation solver.
 */

using namespace heat_equation;

int main()
{
    std::cout << "Example: Gaussian Pulse Diffusion" << std::endl;
    std::cout << "==================================" << std::endl;

    // Set up solver configuration
    SolverConfig config;
    config.alpha = 0.01;       // Thermal diffusivity
    config.dx = 0.01;          // Spatial step
    config.dt = 0.0001;        // Time step
    config.nx = 101;           // Number of grid points
    config.t_final = 0.5;      // Final time
    config.bc = BoundaryCondition::Dirichlet;

    // Create and initialize solver
    HeatSolver1D solver;
    solver.initialize(config);

    // Generate Gaussian initial condition
    // u0(x) = exp(-100 * (x - 0.5)^2)
    std::vector<double> u0 = utils::generateInitialCondition(
        config.nx,
        0.0,    // x_min
        1.0,    // x_max
        [](double x) { return std::exp(-100.0 * (x - 0.5) * (x - 0.5)); }
    );

    // Set initial condition
    solver.setInitialCondition(u0);

    // Set boundary conditions (u = 0 at both ends)
    solver.setBoundaryConditions(0.0, 0.0);

    // Solve the heat equation
    std::cout << "Solving..." << std::endl;
    solver.solve();

    // Export results
    solver.exportSolution("gaussian_diffusion.csv");
    std::cout << "Solution exported to gaussian_diffusion.csv" << std::endl;

    // Print max temperature at final time
    auto solution = solver.getSolution();
    double max_temp = 0.0;
    for (const auto& val : solution) {
        if (val > max_temp) max_temp = val;
    }
    std::cout << "Maximum temperature at t = " << config.t_final << ": " << max_temp << std::endl;

    return 0;
}
