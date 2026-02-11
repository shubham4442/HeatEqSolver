# Heat Equation Solver Documentation

Welcome to the Heat Equation Solver documentation. This solver provides numerical solutions to the heat equation using finite difference methods.

## Documentation Contents

### 1. [Physics and Mathematics](physics_and_math.md)
- Physical background of heat conduction
- Fourier's Law and its derivation
- The heat equation in 1D, 2D, and 3D
- Boundary and initial conditions
- Analytical solutions
- Related equations

### 2. [Numerical Methods](numerical_methods.md)
- Finite difference discretization
- Explicit FTCS scheme
- Von Neumann stability analysis
- **Stability conditions** (critical reading!)
- Convergence theory
- Error sources and verification

### 3. [User Guide](user_guide.md)
- Getting started
- Building the software
- Configuration options
- Setting initial and boundary conditions
- Output and visualization
- Step-by-step examples
- Troubleshooting

### 4. [API Reference](api_reference.md)
- Complete class documentation
- Method signatures and descriptions
- Utility functions
- Code examples

## Quick Reference

### Stability Conditions (IMPORTANT!)

| Dimension | Maximum Time Step |
|-----------|-------------------|
| 1D | Δt ≤ 0.5 × Δx² / α |
| 2D | Δt ≤ 0.25 × min(Δx², Δy²) / α |
| 3D | Δt ≤ (1/6) × min(Δx², Δy², Δz²) / α |

**Violation causes solution explosion!**

### Typical Workflow

```cpp
// 1. Configure
SolverConfig config;
config.alpha = 0.01;
config.dx = 0.01;
config.dt = 0.0001;  // Check stability!
config.nx = 101;
config.t_final = 1.0;
config.bc = BoundaryCondition::Dirichlet;

// 2. Initialize
HeatSolver1D solver;
solver.initialize(config);

// 3. Set conditions
solver.setInitialCondition(u0);
solver.setBoundaryConditions(0.0, 0.0);

// 4. Solve
solver.solve();

// 5. Export
solver.exportSolution("solution.csv");
```

### Visualization

```bash
python scripts/visualize.py solution.csv --dim 1
python scripts/visualize.py solution_2d.csv --dim 2
python scripts/visualize.py solution_3d.csv --dim 3 --slice 10
```

## Document Versions

| Document | Last Updated | Version |
|----------|--------------|---------|
| Physics and Math | 2026-02-11 | 1.0 |
| Numerical Methods | 2026-02-11 | 1.0 |
| User Guide | 2026-02-11 | 1.0 |
| API Reference | 2026-02-11 | 1.0 |

## Additional Resources

- [Heat equation (Wikipedia)](https://en.wikipedia.org/wiki/Heat_equation)
- [Finite difference method (Wikipedia)](https://en.wikipedia.org/wiki/Finite_difference_method)
- LeVeque, R.J. (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations*
