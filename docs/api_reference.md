# API Reference

## Namespace: `heat_equation`

All classes and functions are in the `heat_equation` namespace.

```cpp
using namespace heat_equation;
```

---

## Enumerations

### BoundaryCondition

```cpp
enum class BoundaryCondition {
    Dirichlet,  // Fixed temperature at boundary
    Neumann,    // Fixed heat flux at boundary
    Periodic    // Periodic (wrap-around) boundaries
};
```

---

## Structures

### SolverConfig

Configuration parameters for heat equation solvers.

```cpp
struct SolverConfig {
    double alpha;           // Thermal diffusivity (m²/s)
    double dx;              // Grid spacing in x (m)
    double dy;              // Grid spacing in y (m) - for 2D/3D
    double dz;              // Grid spacing in z (m) - for 3D
    double dt;              // Time step (s)
    size_t nx;              // Number of grid points in x
    size_t ny;              // Number of grid points in y
    size_t nz;              // Number of grid points in z
    double t_final;         // Final simulation time (s)
    BoundaryCondition bc;   // Boundary condition type
    
    // Domain bounds (defaults to [0,1] in each dimension)
    double x_min = 0.0;
    double x_max = 1.0;
    double y_min = 0.0;
    double y_max = 1.0;
    double z_min = 0.0;
    double z_max = 1.0;
};
```

---

## Classes

### Grid1D

1D computational grid.

#### Constructor

```cpp
Grid1D();
Grid1D(size_t n, double x_min, double x_max);
```

#### Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `initialize(size_t n, double x_min, double x_max)` | void | Initialize grid |
| `size()` | size_t | Number of grid points |
| `dx()` | double | Grid spacing |
| `x(size_t i)` | double | Coordinate at index i |
| `coordinates()` | const vector<double>& | All coordinates |

#### Example

```cpp
Grid1D grid(101, 0.0, 1.0);
std::cout << "dx = " << grid.dx() << std::endl;  // 0.01
std::cout << "x[50] = " << grid.x(50) << std::endl;  // 0.5
```

---

### Grid2D

2D computational grid.

#### Constructor

```cpp
Grid2D();
Grid2D(size_t nx, size_t ny, double x_min, double x_max, double y_min, double y_max);
```

#### Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `initialize(...)` | void | Initialize grid |
| `sizeX()`, `sizeY()` | size_t | Grid dimensions |
| `dx()`, `dy()` | double | Grid spacing |
| `x(size_t i)`, `y(size_t j)` | double | Coordinates |

---

### Grid3D

3D computational grid.

#### Constructor

```cpp
Grid3D();
```

#### Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `initialize(size_t nx, size_t ny, size_t nz, ...)` | void | Initialize grid |
| `sizeX()`, `sizeY()`, `sizeZ()` | size_t | Grid dimensions |
| `dx()`, `dy()`, `dz()` | double | Grid spacing |

---

### HeatSolverBase

Abstract base class for all heat solvers.

#### Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `initialize(const SolverConfig&)` | void | Initialize solver |
| `setInitialCondition(const vector<double>&)` | void | Set u(x,0) |
| `setBoundaryConditions(double, double)` | void | Set boundary values |
| `solve()` | void | Run simulation |
| `getSolution()` | vector<double> | Get current solution |
| `exportSolution(const string&)` | void | Export to CSV |
| `getCurrentTime()` | double | Get simulation time |
| `getConfig()` | const SolverConfig& | Get configuration |

---

### HeatSolver1D

1D heat equation solver.

#### Constructor

```cpp
HeatSolver1D();
```

#### Methods

Inherits all methods from `HeatSolverBase`, plus:

| Method | Return Type | Description |
|--------|-------------|-------------|
| `step()` | void | Advance one time step |
| `computeStabilityCondition()` | double | Return r = αΔt/Δx² |
| `getGrid()` | const Grid1D& | Access the grid |

#### Example

```cpp
HeatSolver1D solver;
SolverConfig config;
config.alpha = 0.01;
config.dx = 0.01;
config.dt = 0.0001;
config.nx = 101;
config.t_final = 1.0;
config.bc = BoundaryCondition::Dirichlet;

solver.initialize(config);

std::vector<double> u0(config.nx, 0.0);
u0[50] = 1.0;  // Point source
solver.setInitialCondition(u0);
solver.setBoundaryConditions(0.0, 0.0);

solver.solve();
solver.exportSolution("output.csv");
```

---

### HeatSolver2D

2D heat equation solver.

#### Constructor

```cpp
HeatSolver2D();
```

#### Methods

Inherits all methods from `HeatSolverBase`, plus:

| Method | Return Type | Description |
|--------|-------------|-------------|
| `step()` | void | Advance one time step |
| `computeStabilityCondition()` | double | Return stability parameter |
| `getGrid()` | const Grid2D& | Access the grid |

#### Notes

- Initial condition vector size must be `nx * ny`
- Data is in row-major order: `index = i * ny + j`

---

### HeatSolver3D

3D heat equation solver.

#### Constructor

```cpp
HeatSolver3D();
```

#### Methods

Inherits all methods from `HeatSolverBase`, plus:

| Method | Return Type | Description |
|--------|-------------|-------------|
| `step()` | void | Advance one time step |
| `computeStabilityCondition()` | double | Return stability parameter |
| `getGrid()` | const Grid3D& | Access the grid |
| `exportSlice(const string&, size_t z_index)` | void | Export 2D slice |

#### Notes

- Initial condition vector size must be `nx * ny * nz`
- Data is in row-major order: `index = i * ny * nz + j * nz + k`

---

## Namespace: `heat_equation::utils`

Utility functions.

### File I/O

```cpp
void writeToCsv(const string& filename, 
                const vector<double>& x, 
                const vector<double>& u);

void writeToCsv2D(const string& filename,
                  const vector<vector<double>>& u);

vector<double> readFromFile(const string& filename);
```

### Initial Condition Generation

```cpp
vector<double> generateInitialCondition(
    size_t n,
    double x_min,
    double x_max,
    std::function<double(double)> func
);
```

#### Example

```cpp
auto u0 = utils::generateInitialCondition(
    101, 0.0, 1.0,
    [](double x) { return std::sin(M_PI * x); }
);
```

### Error Computation

```cpp
double computeL2Error(const vector<double>& u, 
                      const vector<double>& u_exact);

double computeMaxError(const vector<double>& u,
                       const vector<double>& u_exact);
```

#### Example

```cpp
auto numerical = solver.getSolution();
auto exact = computeExactSolution(grid, t_final);

double l2_err = utils::computeL2Error(numerical, exact);
double max_err = utils::computeMaxError(numerical, exact);

std::cout << "L2 error: " << l2_err << std::endl;
std::cout << "Max error: " << max_err << std::endl;
```

### Timer

```cpp
class Timer {
public:
    Timer();
    void start();
    void stop();
    double elapsed() const;  // Returns seconds
    void reset();
};
```

#### Example

```cpp
utils::Timer timer;
timer.start();

solver.solve();

timer.stop();
std::cout << "Elapsed: " << timer.elapsed() << " s" << std::endl;
```

---

## Complete Example

```cpp
#include "heat_solver.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>

using namespace heat_equation;

int main() {
    // Configuration
    SolverConfig config;
    config.alpha = 0.01;
    config.dx = 0.01;
    config.dt = 0.0001;
    config.nx = 101;
    config.t_final = 0.5;
    config.bc = BoundaryCondition::Dirichlet;
    config.x_min = 0.0;
    config.x_max = 1.0;

    // Create solver
    HeatSolver1D solver;
    solver.initialize(config);

    // Check stability
    double r = solver.computeStabilityCondition();
    std::cout << "Stability parameter r = " << r;
    std::cout << (r <= 0.5 ? " (stable)" : " (UNSTABLE!)") << std::endl;

    // Initial condition: Gaussian
    const Grid1D& grid = solver.getGrid();
    std::vector<double> u0(config.nx);
    for (size_t i = 0; i < config.nx; ++i) {
        double x = grid.x(i);
        u0[i] = std::exp(-100.0 * (x - 0.5) * (x - 0.5));
    }
    solver.setInitialCondition(u0);

    // Boundary conditions
    solver.setBoundaryConditions(0.0, 0.0);

    // Solve with timing
    utils::Timer timer;
    timer.start();
    solver.solve();
    timer.stop();

    // Output
    solver.exportSolution("solution.csv");
    
    std::cout << "Solved in " << timer.elapsed() << " seconds" << std::endl;
    std::cout << "Final time: " << solver.getCurrentTime() << std::endl;

    return 0;
}
```
