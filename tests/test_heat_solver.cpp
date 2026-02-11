#include "heat_solver.hpp"
#include "grid.hpp"
#include "utils.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

using namespace heat_equation;

/**
 * @brief Test 1D solver initialization
 */
void test_solver1d_init()
{
    HeatSolver1D solver;
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.1;
    config.dt = 0.001;
    config.nx = 11;
    config.t_final = 0.1;
    config.bc = BoundaryCondition::Dirichlet;
    
    solver.initialize(config);
    
    // Should not throw
    std::cout << "test_solver1d_init: PASSED" << std::endl;
}

/**
 * @brief Test stability condition calculation
 */
void test_stability_condition()
{
    HeatSolver1D solver;
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.1;
    config.dt = 0.001;
    config.nx = 11;
    config.t_final = 0.1;
    config.bc = BoundaryCondition::Dirichlet;
    
    solver.initialize(config);
    
    // r = alpha * dt / dx^2 = 0.01 * 0.001 / 0.01 = 0.001
    double r = solver.computeStabilityCondition();
    assert(r < 0.5);  // Should be stable
    
    std::cout << "test_stability_condition: PASSED (r = " << r << ")" << std::endl;
}

/**
 * @brief Test initial condition setting
 */
void test_initial_condition()
{
    HeatSolver1D solver;
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.1;
    config.dt = 0.001;
    config.nx = 11;
    config.t_final = 0.1;
    config.bc = BoundaryCondition::Dirichlet;
    
    solver.initialize(config);
    
    std::vector<double> u0(config.nx, 1.0);
    solver.setInitialCondition(u0);
    
    auto solution = solver.getSolution();
    assert(solution.size() == config.nx);
    assert(std::abs(solution[5] - 1.0) < 1e-10);
    
    std::cout << "test_initial_condition: PASSED" << std::endl;
}

/**
 * @brief Test single time step
 */
void test_single_step()
{
    HeatSolver1D solver;
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.1;
    config.dt = 0.001;
    config.nx = 11;
    config.t_final = 0.001;  // Single step
    config.bc = BoundaryCondition::Dirichlet;
    
    solver.initialize(config);
    
    // Set a simple initial condition
    std::vector<double> u0(config.nx, 0.0);
    u0[5] = 1.0;  // Delta function in middle
    solver.setInitialCondition(u0);
    solver.setBoundaryConditions(0.0, 0.0);
    
    solver.step();
    auto solution = solver.getSolution();
    
    // After one step, heat should spread to neighbors
    assert(solution[5] < 1.0);  // Center should decrease
    assert(solution[4] > 0.0);  // Neighbors should increase
    assert(solution[6] > 0.0);
    
    std::cout << "test_single_step: PASSED" << std::endl;
}

/**
 * @brief Test boundary conditions
 */
void test_boundary_conditions()
{
    HeatSolver1D solver;
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.1;
    config.dt = 0.001;
    config.nx = 11;
    config.t_final = 0.01;
    config.bc = BoundaryCondition::Dirichlet;
    
    solver.initialize(config);
    
    std::vector<double> u0(config.nx, 1.0);
    solver.setInitialCondition(u0);
    solver.setBoundaryConditions(0.0, 0.0);
    
    solver.solve();
    auto solution = solver.getSolution();
    
    // Dirichlet BC: boundaries should be 0
    assert(std::abs(solution[0]) < 1e-10);
    assert(std::abs(solution[config.nx - 1]) < 1e-10);
    
    std::cout << "test_boundary_conditions: PASSED" << std::endl;
}

/**
 * @brief Test energy conservation (for periodic BC)
 */
void test_energy_conservation()
{
    HeatSolver1D solver;
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.1;
    config.dt = 0.0001;  // Very small dt for stability
    config.nx = 21;
    config.t_final = 0.01;
    config.bc = BoundaryCondition::Periodic;
    
    solver.initialize(config);
    
    std::vector<double> u0(config.nx);
    for (size_t i = 0; i < config.nx; ++i) {
        u0[i] = std::sin(2.0 * M_PI * i / (config.nx - 1));
    }
    solver.setInitialCondition(u0);
    
    // Calculate initial total heat
    double initial_sum = 0.0;
    for (const auto& val : u0) initial_sum += val;
    
    solver.solve();
    auto solution = solver.getSolution();
    
    // Calculate final total heat
    double final_sum = 0.0;
    for (const auto& val : solution) final_sum += val;
    
    // For periodic BC, total heat should be approximately conserved
    assert(std::abs(final_sum - initial_sum) < 0.1);
    
    std::cout << "test_energy_conservation: PASSED" << std::endl;
}

/**
 * @brief Test 2D solver
 */
void test_solver2d()
{
    HeatSolver2D solver;
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.1;
    config.dt = 0.0001;
    config.nx = 11;
    config.ny = 11;
    config.t_final = 0.001;
    config.bc = BoundaryCondition::Dirichlet;
    
    solver.initialize(config);
    
    std::vector<double> u0(config.nx * config.ny, 0.0);
    u0[5 * config.ny + 5] = 1.0;  // Center point
    solver.setInitialCondition(u0);
    solver.setBoundaryConditions(0.0, 0.0);
    
    solver.solve();
    auto solution = solver.getSolution();
    
    assert(solution.size() == config.nx * config.ny);
    
    std::cout << "test_solver2d: PASSED" << std::endl;
}

int main()
{
    std::cout << "Running Heat Equation Solver Tests" << std::endl;
    std::cout << "===================================" << std::endl;

    test_solver1d_init();
    test_stability_condition();
    test_initial_condition();
    test_single_step();
    test_boundary_conditions();
    test_energy_conservation();
    test_solver2d();

    std::cout << "===================================" << std::endl;
    std::cout << "All tests passed!" << std::endl;

    return 0;
}
