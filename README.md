# Heat Equation Solver

A C++ library for solving the heat equation using finite difference methods, with Python visualization tools.

## Overview

This project provides solvers for the heat equation (also known as the diffusion equation):

```
∂u/∂t = α ∇²u
```

where:
- `u` is the temperature distribution
- `α` is the thermal diffusivity
- `∇²` is the Laplacian operator

## Features

- **1D Heat Solver**: Explicit finite difference method for 1D problems
- **2D Heat Solver**: Extension to 2D problems with proper stability handling
- **3D Heat Solver**: Full 3D solver with slice export capability
- **Grid Classes**: Flexible grid generation for 1D, 2D, and 3D domains
- **Utility Functions**: File I/O, error computation, timing
- **Python Visualization**: Scripts for 1D, 2D, and 3D solution visualization

## Project Structure

```
heat_equation_solver/
├── include/
│   ├── heat_solver.hpp    # Main solver classes (1D, 2D, 3D)
│   ├── grid.hpp           # Grid classes (Grid1D, Grid2D, Grid3D)
│   └── utils.hpp          # Utility functions
├── src/
│   ├── heat_solver.cpp    # Solver implementations
│   ├── grid.cpp           # Grid implementations
│   ├── utils.cpp          # Utility implementations
│   └── main.cpp           # Main program
├── tests/
│   ├── test_heat_solver.cpp
│   └── test_grid.cpp
├── examples/
│   ├── gaussian_diffusion.cpp  # 1D Gaussian pulse example
│   └── heat_3d.cpp             # 3D heat diffusion example
├── scripts/
│   ├── visualize.py           # Visualization script
│   ├── generate_test_data.py  # Generate synthetic test data
│   └── requirements.txt       # Python dependencies
├── docs/
├── CMakeLists.txt
├── Makefile
└── README.md
```

## Building

### Requirements

- C++17 compatible compiler (g++, clang++)
- Make or CMake 3.14+

### Build with Make

```bash
make all       # Build everything
make test      # Run tests
make run       # Run the main 1D solver
make run3d     # Run the 3D example
make clean     # Clean build files
make help      # Show available targets
```

### Build with CMake

```bash
mkdir build && cd build
cmake ..
make
ctest          # Run tests
```

## Usage

### C++ API

```cpp
#include "heat_solver.hpp"

using namespace heat_equation;

int main() {
    // Configure solver
    SolverConfig config;
    config.alpha = 0.01;        // Thermal diffusivity
    config.dx = 0.01;           // Spatial step
    config.dt = 0.0001;         // Time step
    config.nx = 101;            // Grid points
    config.t_final = 1.0;       // Final time
    config.bc = BoundaryCondition::Dirichlet;
    config.x_min = 0.0;
    config.x_max = 1.0;

    // Create and run 1D solver
    HeatSolver1D solver;
    solver.initialize(config);
    
    // Use the grid for initial condition
    const Grid1D& grid = solver.getGrid();
    std::vector<double> u0(config.nx);
    for (size_t i = 0; i < config.nx; ++i) {
        u0[i] = std::exp(-100.0 * pow(grid.x(i) - 0.5, 2));
    }
    
    solver.setInitialCondition(u0);
    solver.setBoundaryConditions(0.0, 0.0);
    solver.solve();
    solver.exportSolution("solution.csv");

    return 0;
}
```

### 3D Example

```cpp
// Configure for 3D
SolverConfig config;
config.alpha = 0.01;
config.nx = 21; config.ny = 21; config.nz = 21;
config.dt = 0.00001;  // Smaller dt for 3D stability
config.t_final = 0.1;

HeatSolver3D solver;
solver.initialize(config);
// ... set initial conditions ...
solver.solve();
solver.exportSolution("solution_3d.csv");
solver.exportSlice("slice.csv", config.nz / 2);  // Export middle slice
```

## Python Visualization

### Installation

```bash
pip install -r scripts/requirements.txt
```

### Visualizing Solutions

```bash
# 1D solution
python scripts/visualize.py solution.csv --dim 1

# 2D solution (contour + surface plot)
python scripts/visualize.py solution_2d.csv --dim 2

# 2D solution (heatmap only)
python scripts/visualize.py solution_2d.csv --dim 2 --type heatmap

# 3D solution (slice at z=10)
python scripts/visualize.py solution_3d.csv --dim 3 --slice 10

# 3D solution (isosurface at 50% of max temperature)
python scripts/visualize.py solution_3d.csv --dim 3 --iso 0.5

# Save plot to file
python scripts/visualize.py solution.csv --dim 1 --save output.png
```

### Generate Test Data

```bash
python scripts/generate_test_data.py
```

## Numerical Method

### Explicit Finite Difference Scheme

**1D:**
```
u[i]^{n+1} = u[i]^n + r * (u[i+1]^n - 2*u[i]^n + u[i-1]^n)
```

**2D:**
```
u[i,j]^{n+1} = u[i,j]^n + rx*(u[i+1,j] - 2*u[i,j] + u[i-1,j])
                        + ry*(u[i,j+1] - 2*u[i,j] + u[i,j-1])
```

**3D:**
```
u[i,j,k]^{n+1} = u[i,j,k]^n + rx*(u[i+1,j,k] - 2*u + u[i-1,j,k])
                            + ry*(u[i,j+1,k] - 2*u + u[i,j-1,k])
                            + rz*(u[i,j,k+1] - 2*u + u[i,j,k-1])
```

where `r = α * dt / dx²`.

### Stability Conditions

| Dimension | Stability Requirement |
|-----------|----------------------|
| 1D        | r ≤ 0.5              |
| 2D        | r ≤ 0.25             |
| 3D        | r ≤ 1/6 ≈ 0.167      |

## Boundary Conditions

- **Dirichlet**: Fixed temperature at boundaries
- **Neumann**: Fixed heat flux at boundaries  
- **Periodic**: Periodic boundary conditions

## Grid Classes

The `Grid1D`, `Grid2D`, and `Grid3D` classes provide:
- Uniform grid generation
- Coordinate access via `x(i)`, `y(j)`, `z(k)`
- Grid spacing via `dx()`, `dy()`, `dz()`
- Coordinate vectors via `coordinates()`

## License

MIT License

## TODO

- [ ] Implement implicit methods (Crank-Nicolson)
- [ ] Add parallel computation support (OpenMP)
- [ ] Add VTK output for 3D visualization in ParaView
- [ ] Implement adaptive time stepping
