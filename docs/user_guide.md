# User Guide: Heat Equation Solver

## Table of Contents

1. [Getting Started](#getting-started)
2. [Building the Software](#building-the-software)
3. [Running Your First Simulation](#running-your-first-simulation)
4. [Configuration Options](#configuration-options)
5. [Setting Initial Conditions](#setting-initial-conditions)
6. [Boundary Conditions](#boundary-conditions)
7. [Output and Visualization](#output-and-visualization)
8. [Examples](#examples)
9. [Troubleshooting](#troubleshooting)

## Getting Started

### Prerequisites

**C++ Solver:**
- C++17 compatible compiler (g++ 7+, clang++ 5+)
- Make or CMake 3.14+

**Python Visualization:**
- Python 3.7+
- NumPy, Matplotlib

### Quick Start

```bash
# Clone/navigate to the project
cd heat_equation_solver

# Build
make all

# Run 1D example
./build/heat_solver

# Visualize results
pip install -r scripts/requirements.txt
python scripts/visualize.py solution.csv --dim 1
```

## Building the Software

### Using Make (Recommended)

```bash
# Build all targets
make all

# Build specific targets
make build/heat_solver        # Main solver
make build/test_heat_solver   # Tests
make build/example_3d         # 3D example

# Show available targets
make help
```

### Using CMake

```bash
mkdir build && cd build
cmake ..
make -j4

# Run tests
ctest
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `CXX` | g++ | C++ compiler |
| `CXXFLAGS` | -std=c++17 -O2 | Compiler flags |
| `-DCMAKE_BUILD_TYPE` | Release | CMake build type |

## Running Your First Simulation

### Step 1: Create Configuration

```cpp
#include "heat_solver.hpp"

using namespace heat_equation;

SolverConfig config;
config.alpha = 0.01;      // Thermal diffusivity (m²/s)
config.dx = 0.01;         // Grid spacing (m)
config.dt = 0.0001;       // Time step (s)
config.nx = 101;          // Number of grid points
config.t_final = 1.0;     // Simulation duration (s)
config.bc = BoundaryCondition::Dirichlet;

// Domain bounds
config.x_min = 0.0;
config.x_max = 1.0;
```

### Step 2: Initialize Solver

```cpp
HeatSolver1D solver;
solver.initialize(config);
```

### Step 3: Set Initial Condition

```cpp
// Get grid for coordinate access
const Grid1D& grid = solver.getGrid();

// Create initial condition
std::vector<double> u0(config.nx);
for (size_t i = 0; i < config.nx; ++i) {
    double x = grid.x(i);
    u0[i] = std::exp(-100.0 * (x - 0.5) * (x - 0.5));  // Gaussian
}

solver.setInitialCondition(u0);
```

### Step 4: Set Boundary Conditions

```cpp
solver.setBoundaryConditions(0.0, 0.0);  // u=0 at both ends
```

### Step 5: Solve and Export

```cpp
solver.solve();
solver.exportSolution("solution.csv");
```

## Configuration Options

### SolverConfig Structure

| Field | Type | Description | Units |
|-------|------|-------------|-------|
| `alpha` | double | Thermal diffusivity | m²/s |
| `dx` | double | Grid spacing in x | m |
| `dy` | double | Grid spacing in y (2D/3D) | m |
| `dz` | double | Grid spacing in z (3D) | m |
| `dt` | double | Time step | s |
| `nx` | size_t | Grid points in x | - |
| `ny` | size_t | Grid points in y | - |
| `nz` | size_t | Grid points in z | - |
| `t_final` | double | Simulation end time | s |
| `bc` | BoundaryCondition | Boundary type | - |
| `x_min`, `x_max` | double | Domain bounds | m |
| `y_min`, `y_max` | double | Domain bounds (2D/3D) | m |
| `z_min`, `z_max` | double | Domain bounds (3D) | m |

### Choosing Parameters

#### Grid Spacing (dx)

```cpp
// Rule of thumb: 5-10 points per feature
double feature_width = 0.1;  // Width of initial pulse
config.dx = feature_width / 10;
config.nx = static_cast<size_t>((config.x_max - config.x_min) / config.dx) + 1;
```

#### Time Step (dt)

The time step must satisfy the stability condition:

```cpp
// 1D: dt <= 0.5 * dx² / alpha
double dt_max_1d = 0.5 * config.dx * config.dx / config.alpha;

// 2D: dt <= 0.25 * min(dx², dy²) / alpha
double dt_max_2d = 0.25 * std::min(dx*dx, dy*dy) / config.alpha;

// 3D: dt <= (1/6) * min(dx², dy², dz²) / alpha
double dt_max_3d = (1.0/6.0) * std::min({dx*dx, dy*dy, dz*dz}) / config.alpha;

// Use 90% of maximum for safety
config.dt = 0.9 * dt_max;
```

#### Thermal Diffusivity (alpha)

Common values:

| Material | α (m²/s) |
|----------|----------|
| Copper | 1.1 × 10⁻⁴ |
| Steel | 1.2 × 10⁻⁵ |
| Glass | 3.4 × 10⁻⁷ |
| Water | 1.4 × 10⁻⁷ |

For dimensionless simulations, use α = 1.

## Setting Initial Conditions

### Using Lambda Functions

```cpp
// Gaussian pulse
auto gaussian = [](double x) {
    return std::exp(-100.0 * (x - 0.5) * (x - 0.5));
};

// Sinusoidal
auto sine = [](double x) {
    return std::sin(M_PI * x);
};

// Step function
auto step = [](double x) {
    return (x < 0.5) ? 1.0 : 0.0;
};

// Generate initial condition
std::vector<double> u0 = utils::generateInitialCondition(
    config.nx, config.x_min, config.x_max, gaussian
);
```

### Manual Generation

```cpp
std::vector<double> u0(config.nx);
const Grid1D& grid = solver.getGrid();

for (size_t i = 0; i < config.nx; ++i) {
    double x = grid.x(i);
    // Your custom function
    u0[i] = ...;
}
```

### 2D Initial Conditions

```cpp
std::vector<double> u0(config.nx * config.ny);
const Grid2D& grid = solver.getGrid();

for (size_t i = 0; i < config.nx; ++i) {
    for (size_t j = 0; j < config.ny; ++j) {
        double x = grid.x(i);
        double y = grid.y(j);
        size_t idx = i * config.ny + j;  // Row-major order
        
        // 2D Gaussian
        double r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);
        u0[idx] = std::exp(-50.0 * r2);
    }
}
```

## Boundary Conditions

### Available Types

| Type | Description | Use Case |
|------|-------------|----------|
| `Dirichlet` | Fixed temperature | Heat reservoir |
| `Neumann` | Fixed flux | Insulated wall |
| `Periodic` | Wrapping | Ring geometry |

### Setting Boundary Conditions

```cpp
// Dirichlet: u(0) = 1.0, u(L) = 0.0
config.bc = BoundaryCondition::Dirichlet;
solver.setBoundaryConditions(1.0, 0.0);

// Neumann: du/dx(0) = 0, du/dx(L) = 0 (insulated)
config.bc = BoundaryCondition::Neumann;
solver.setBoundaryConditions(0.0, 0.0);

// Periodic: u(0) = u(L)
config.bc = BoundaryCondition::Periodic;
// No need to call setBoundaryConditions for periodic
```

## Output and Visualization

### Export Formats

**1D Output (CSV):**
```
x,u
0.0,0.0
0.01,0.00135
...
```

**2D Output (CSV):**
```
x,y,u
0.0,0.0,0.0
0.0,0.1,0.05
...
```

**3D Output (CSV):**
```
x,y,z,u
0.0,0.0,0.0,0.0
...
```

### Python Visualization

```bash
# Install dependencies
pip install -r scripts/requirements.txt

# 1D: Line plot
python scripts/visualize.py solution.csv --dim 1

# 2D: Contour + surface
python scripts/visualize.py solution_2d.csv --dim 2

# 2D: Heatmap only
python scripts/visualize.py solution_2d.csv --dim 2 --type heatmap

# 3D: Slice at z-index 10
python scripts/visualize.py solution_3d.csv --dim 3 --slice 10

# 3D: Slice along y-axis
python scripts/visualize.py solution_3d.csv --dim 3 --slice 5 --slice-axis y

# 3D: Isosurface at 50% of max temperature
python scripts/visualize.py solution_3d.csv --dim 3 --iso 0.5

# Save to file
python scripts/visualize.py solution.csv --dim 1 --save result.png

# Custom title
python scripts/visualize.py solution.csv --dim 1 --title "My Simulation"
```

### Programmatic Access

```cpp
// Get solution as vector
std::vector<double> solution = solver.getSolution();

// Access specific values
double center_temp = solution[config.nx / 2];

// For 3D: export just a slice
solver.exportSlice("slice.csv", config.nz / 2);
```

## Examples

### Example 1: Cooling of a Hot Rod

A rod initially at 100°C, cooled at both ends to 0°C.

```cpp
SolverConfig config;
config.alpha = 1e-5;     // Steel
config.dx = 0.001;       // 1 mm resolution
config.dt = 0.001;       // Calculated for stability
config.nx = 101;         // 10 cm rod
config.t_final = 60.0;   // 1 minute
config.bc = BoundaryCondition::Dirichlet;

HeatSolver1D solver;
solver.initialize(config);

// Initially at 100°C
std::vector<double> u0(config.nx, 100.0);
solver.setInitialCondition(u0);
solver.setBoundaryConditions(0.0, 0.0);  // Cooled ends

solver.solve();
solver.exportSolution("cooling_rod.csv");
```

### Example 2: Heat Spreading on a Plate

2D heat diffusion from a central hot spot.

```cpp
SolverConfig config;
config.alpha = 1e-4;     // Copper plate
config.nx = 51; config.ny = 51;
config.dt = 0.00001;
config.t_final = 0.5;
config.bc = BoundaryCondition::Dirichlet;

HeatSolver2D solver;
solver.initialize(config);

// Hot spot in center
std::vector<double> u0(config.nx * config.ny, 0.0);
for (size_t i = 20; i < 30; ++i) {
    for (size_t j = 20; j < 30; ++j) {
        u0[i * config.ny + j] = 100.0;
    }
}

solver.setInitialCondition(u0);
solver.setBoundaryConditions(0.0, 0.0);
solver.solve();
solver.exportSolution("plate_heat.csv");
```

### Example 3: 3D Heat Diffusion in a Cube

```cpp
SolverConfig config;
config.alpha = 1e-4;
config.nx = 21; config.ny = 21; config.nz = 21;
config.dt = 1e-6;
config.t_final = 0.1;
config.bc = BoundaryCondition::Dirichlet;

HeatSolver3D solver;
solver.initialize(config);

// Spherical hot region
std::vector<double> u0(config.nx * config.ny * config.nz, 0.0);
for (size_t i = 0; i < config.nx; ++i) {
    for (size_t j = 0; j < config.ny; ++j) {
        for (size_t k = 0; k < config.nz; ++k) {
            double x = i * config.dx - 0.5;
            double y = j * config.dx - 0.5;
            double z = k * config.dx - 0.5;
            double r2 = x*x + y*y + z*z;
            if (r2 < 0.01) {
                size_t idx = i * config.ny * config.nz + j * config.nz + k;
                u0[idx] = 100.0;
            }
        }
    }
}

solver.setInitialCondition(u0);
solver.setBoundaryConditions(0.0, 0.0);
solver.solve();
solver.exportSolution("cube_heat.csv");
solver.exportSlice("cube_slice.csv", config.nz / 2);
```

## Troubleshooting

### Solution Explodes (Goes to Infinity)

**Cause**: Stability condition violated (dt too large)

**Solution**: Reduce dt or check stability:
```cpp
double r = config.alpha * config.dt / (config.dx * config.dx);
std::cout << "r = " << r << " (must be <= 0.5 for 1D)" << std::endl;
```

### Solution Doesn't Change

**Cause**: t_final too small or dt too large

**Solution**: 
```cpp
size_t num_steps = config.t_final / config.dt;
std::cout << "Number of time steps: " << num_steps << std::endl;
// Should be > 0
```

### "Initial condition size does not match"

**Cause**: Wrong vector size for initial condition

**Solution**:
```cpp
// 1D: size = nx
// 2D: size = nx * ny
// 3D: size = nx * ny * nz
```

### Poor Resolution

**Cause**: Grid too coarse

**Solution**: Increase nx (and ny, nz for multi-D)

### Simulation Too Slow

**Causes and Solutions**:

1. **Grid too fine**: Use coarser grid if accuracy allows
2. **Too many time steps**: Increase dt (within stability limit)
3. **3D simulation**: Reduce grid size or use slices

### Visualization Errors

**"File not found"**: Check file path and working directory

**"ValueError: could not broadcast"**: Check CSV format matches expected dimension

## Performance Tips

1. **Choose dt at stability limit**: Maximizes efficiency while stable
2. **Use appropriate precision**: double is usually sufficient
3. **Profile your code**: Use `utils::Timer` to measure performance
4. **Consider parallel computing**: Future feature with OpenMP

## Getting Help

1. Check the documentation in `docs/`
2. Run tests: `make test`
3. Generate test data: `python scripts/generate_test_data.py`
4. Review examples in `examples/`
