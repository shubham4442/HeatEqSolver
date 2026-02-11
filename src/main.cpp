#include "heat_solver.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>
#include <vector>

using namespace heat_equation;

int main(int argc, char* argv[])
{
    std::cout << "Heat Equation Solver" << std::endl;
    std::cout << "====================" << std::endl;

    // Solver configuration
    SolverConfig config;
    config.alpha = 0.01;        // Thermal diffusivity
    config.dx = 0.01;           // Spatial step
    config.dt = 0.0001;         // Time step (stable: r = alpha*dt/dx^2 = 0.1 <= 0.5)
    config.nx = 101;            // Number of grid points
    config.t_final = 1.0;       // Final time
    config.bc = BoundaryCondition::Dirichlet;
    
    // Domain bounds
    config.x_min = 0.0;
    config.x_max = 1.0;

    // Create solver instance
    HeatSolver1D solver;

    // Initialize solver (this also creates the grid)
    solver.initialize(config);

    // Generate initial condition using the grid coordinates
    const Grid1D& grid = solver.getGrid();
    std::vector<double> u0(config.nx);
    for (size_t i = 0; i < config.nx; ++i) {
        double x = grid.x(i);
        // Gaussian pulse centered at x = 0.5
        u0[i] = std::exp(-100.0 * (x - 0.5) * (x - 0.5));
    }

    // Set initial condition
    solver.setInitialCondition(u0);

    // Set boundary conditions (Dirichlet: u = 0 at boundaries)
    solver.setBoundaryConditions(0.0, 0.0);

    // Time the solve
    utils::Timer timer;
    timer.start();

    // Solve the heat equation
    solver.solve();

    timer.stop();

    // Get solution
    auto solution = solver.getSolution();

    // Export solution
    solver.exportSolution("solution.csv");

    // Print summary
    std::cout << "Solver completed successfully!" << std::endl;
    std::cout << "  Grid points: " << config.nx << std::endl;
    std::cout << "  Domain:      [" << config.x_min << ", " << config.x_max << "]" << std::endl;
    std::cout << "  Grid spacing: " << grid.dx() << std::endl;
    std::cout << "  Final time:  " << config.t_final << std::endl;
    std::cout << "  Time steps:  " << static_cast<size_t>(config.t_final / config.dt) << std::endl;
    std::cout << "  Wall time:   " << timer.elapsed() << " seconds" << std::endl;
    std::cout << "  Solution exported to: solution.csv" << std::endl;

    // Print solution at a few points
    std::cout << "\nSolution at selected points:" << std::endl;
    std::cout << "  u(" << grid.x(0) << ") = " << solution[0] << std::endl;
    std::cout << "  u(" << grid.x(config.nx / 2) << ") = " << solution[config.nx / 2] << std::endl;
    std::cout << "  u(" << grid.x(config.nx - 1) << ") = " << solution[config.nx - 1] << std::endl;

    std::cout << "\nVisualize with:" << std::endl;
    std::cout << "  python scripts/visualize.py solution.csv --dim 1" << std::endl;

    return 0;
}
