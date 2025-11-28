# C++ Heat Equation Solver Project

A C++ based learning project to implement Numeric methods for Solving PDEs

## Features (Planned)
### Core Solvers
- 1D Heat equation: Explicit (FTCS)
- 1D Heat equation: Implicit (BTCS)
- 2D Heat equation: Jacobi, Gauss-Seidel, SOR
- Crank–Nicolson schemes
- Sparse matrix solvers and iterative methods
### Tools & Enhancements
- Multithreaded solvers (std::thread / std::async / thread pool)
- Benchmarking with Google Benchmark
- Unit testing with GoogleTest
- Python visualization using NumPy/Matplotlib
- Config-driven simulations (JSON or YAML)
### Future includes
- FE simulations 

## Build instructions
```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build 
