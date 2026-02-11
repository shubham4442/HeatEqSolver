#include "heat_solver.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>
#include <vector>

/**
 * @brief Example: 3D Heat diffusion with a hot spot in the center
 */

using namespace heat_equation;

int main()
{
    std::cout << "Example: 3D Heat Diffusion" << std::endl;
    std::cout << "==========================" << std::endl;

    // Set up solver configuration for 3D
    SolverConfig config;
    config.alpha = 0.01;        // Thermal diffusivity
    config.nx = 21;             // Grid points in x
    config.ny = 21;             // Grid points in y
    config.nz = 21;             // Grid points in z
    config.dt = 0.00001;        // Small time step for 3D stability
    config.t_final = 0.1;       // Final time
    config.bc = BoundaryCondition::Dirichlet;
    
    // Domain bounds
    config.x_min = 0.0; config.x_max = 1.0;
    config.y_min = 0.0; config.y_max = 1.0;
    config.z_min = 0.0; config.z_max = 1.0;
    
    // Calculate dx, dy, dz
    config.dx = (config.x_max - config.x_min) / (config.nx - 1);
    config.dy = (config.y_max - config.y_min) / (config.ny - 1);
    config.dz = (config.z_max - config.z_min) / (config.nz - 1);

    // Create and initialize solver
    HeatSolver3D solver;
    solver.initialize(config);

    // Generate 3D initial condition: Gaussian hot spot at center
    std::vector<double> u0(config.nx * config.ny * config.nz);
    const Grid3D& grid = solver.getGrid();
    
    double sigma = 0.15;
    for (size_t i = 0; i < config.nx; ++i) {
        double x = config.x_min + i * grid.dx();
        for (size_t j = 0; j < config.ny; ++j) {
            double y = config.y_min + j * grid.dy();
            for (size_t k = 0; k < config.nz; ++k) {
                double z = config.z_min + k * grid.dz();
                
                double r2 = (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) + (z - 0.5)*(z - 0.5);
                size_t idx = i * config.ny * config.nz + j * config.nz + k;
                u0[idx] = std::exp(-r2 / (sigma * sigma));
            }
        }
    }

    // Set initial condition
    solver.setInitialCondition(u0);

    // Set boundary conditions (u = 0 on all faces)
    solver.setBoundaryConditions(0.0, 0.0);

    // Solve the heat equation
    std::cout << "Solving 3D heat equation..." << std::endl;
    std::cout << "  Grid: " << config.nx << " x " << config.ny << " x " << config.nz << std::endl;
    std::cout << "  Total points: " << config.nx * config.ny * config.nz << std::endl;
    
    utils::Timer timer;
    timer.start();
    
    solver.solve();
    
    timer.stop();

    // Export full 3D solution
    solver.exportSolution("solution_3d.csv");
    
    // Export middle slice for easier visualization
    solver.exportSlice("solution_3d_slice.csv", config.nz / 2);
    
    std::cout << "Solution computed in " << timer.elapsed() << " seconds" << std::endl;
    std::cout << "Exported to: solution_3d.csv" << std::endl;
    std::cout << "Middle slice: solution_3d_slice.csv" << std::endl;

    // Find max temperature
    auto solution = solver.getSolution();
    double max_temp = 0.0;
    for (const auto& val : solution) {
        if (val > max_temp) max_temp = val;
    }
    std::cout << "Maximum temperature at t = " << config.t_final << ": " << max_temp << std::endl;

    std::cout << "\nVisualize with:" << std::endl;
    std::cout << "  python scripts/visualize.py solution_3d.csv --dim 3 --slice 10" << std::endl;
    std::cout << "  python scripts/visualize.py solution_3d_slice.csv --dim 2" << std::endl;

    return 0;
}
