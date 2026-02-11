#include "grid.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace heat_equation;

/**
 * @brief Test 1D grid construction
 */
void test_grid1d_construction()
{
    Grid1D grid(11, 0.0, 1.0);
    
    assert(grid.size() == 11);
    assert(std::abs(grid.x(0) - 0.0) < 1e-10);
    assert(std::abs(grid.x(10) - 1.0) < 1e-10);
    
    std::cout << "test_grid1d_construction: PASSED" << std::endl;
}

/**
 * @brief Test 1D grid spacing calculation
 */
void test_grid1d_spacing()
{
    Grid1D grid(11, 0.0, 1.0);
    
    double expected_dx = 0.1;
    assert(std::abs(grid.dx() - expected_dx) < 1e-10);
    
    // Check that spacing is uniform
    const auto& coords = grid.coordinates();
    for (size_t i = 1; i < coords.size(); ++i) {
        double actual_dx = coords[i] - coords[i-1];
        assert(std::abs(actual_dx - expected_dx) < 1e-10);
    }
    
    std::cout << "test_grid1d_spacing: PASSED" << std::endl;
}

/**
 * @brief Test 2D grid construction
 */
void test_grid2d_construction()
{
    Grid2D grid(11, 21, 0.0, 1.0, 0.0, 2.0);
    
    assert(grid.sizeX() == 11);
    assert(grid.sizeY() == 21);
    assert(std::abs(grid.x(0) - 0.0) < 1e-10);
    assert(std::abs(grid.x(10) - 1.0) < 1e-10);
    assert(std::abs(grid.y(0) - 0.0) < 1e-10);
    assert(std::abs(grid.y(20) - 2.0) < 1e-10);
    assert(std::abs(grid.dx() - 0.1) < 1e-10);
    assert(std::abs(grid.dy() - 0.1) < 1e-10);
    
    std::cout << "test_grid2d_construction: PASSED" << std::endl;
}

/**
 * @brief Test 3D grid construction
 */
void test_grid3d_construction()
{
    Grid3D grid;
    grid.initialize(11, 11, 11, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    
    assert(grid.sizeX() == 11);
    assert(grid.sizeY() == 11);
    assert(grid.sizeZ() == 11);
    assert(std::abs(grid.dx() - 0.1) < 1e-10);
    assert(std::abs(grid.dy() - 0.1) < 1e-10);
    assert(std::abs(grid.dz() - 0.1) < 1e-10);
    
    std::cout << "test_grid3d_construction: PASSED" << std::endl;
}

/**
 * @brief Test default grid initialization
 */
void test_grid1d_default()
{
    Grid1D grid;
    assert(grid.size() == 0);
    assert(grid.dx() == 0.0);
    
    // Re-initialize
    grid.initialize(21, -1.0, 1.0);
    assert(grid.size() == 21);
    assert(std::abs(grid.dx() - 0.1) < 1e-10);
    assert(std::abs(grid.x(0) - (-1.0)) < 1e-10);
    assert(std::abs(grid.x(20) - 1.0) < 1e-10);
    
    std::cout << "test_grid1d_default: PASSED" << std::endl;
}

int main()
{
    std::cout << "Running Grid Tests" << std::endl;
    std::cout << "==================" << std::endl;

    test_grid1d_construction();
    test_grid1d_spacing();
    test_grid2d_construction();
    test_grid3d_construction();
    test_grid1d_default();

    std::cout << "==================" << std::endl;
    std::cout << "All tests passed!" << std::endl;

    return 0;
}
