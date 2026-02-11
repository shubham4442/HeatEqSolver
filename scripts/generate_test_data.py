#!/usr/bin/env python3
"""
Generate test data and demonstrate visualization capabilities.

This script creates synthetic heat equation solutions for testing the
visualization tools without needing to run the C++ solver.
"""

import numpy as np
import os


def generate_1d_solution(filename, nx=101, t=0.5, alpha=0.01):
    """Generate 1D Gaussian diffusion solution."""
    x = np.linspace(0, 1, nx)
    
    # Analytical solution for Gaussian initial condition
    # u(x,t) = 1/sqrt(1 + 4*alpha*t) * exp(-(x-0.5)^2 / (sigma^2 + 4*alpha*t))
    sigma = 0.1
    denom = sigma**2 + 4 * alpha * t
    u = (sigma / np.sqrt(denom)) * np.exp(-(x - 0.5)**2 / denom)
    
    with open(filename, 'w') as f:
        f.write("x,u\n")
        for xi, ui in zip(x, u):
            f.write(f"{xi},{ui}\n")
    
    print(f"Generated 1D solution: {filename}")


def generate_2d_solution(filename, nx=51, ny=51, t=0.5, alpha=0.01):
    """Generate 2D Gaussian diffusion solution."""
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # 2D Gaussian diffusion
    sigma = 0.1
    denom = sigma**2 + 4 * alpha * t
    r2 = (X - 0.5)**2 + (Y - 0.5)**2
    U = (sigma**2 / denom) * np.exp(-r2 / denom)
    
    with open(filename, 'w') as f:
        f.write("x,y,u\n")
        for i in range(nx):
            for j in range(ny):
                f.write(f"{X[i,j]},{Y[i,j]},{U[i,j]}\n")
    
    print(f"Generated 2D solution: {filename}")


def generate_3d_solution(filename, nx=21, ny=21, nz=21, t=0.5, alpha=0.01):
    """Generate 3D Gaussian diffusion solution."""
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    z = np.linspace(0, 1, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # 3D Gaussian diffusion
    sigma = 0.15
    denom = sigma**2 + 4 * alpha * t
    r2 = (X - 0.5)**2 + (Y - 0.5)**2 + (Z - 0.5)**2
    U = (sigma**3 / denom**1.5) * np.exp(-r2 / denom)
    
    with open(filename, 'w') as f:
        f.write("x,y,z,u\n")
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    f.write(f"{X[i,j,k]},{Y[i,j,k]},{Z[i,j,k]},{U[i,j,k]}\n")
    
    print(f"Generated 3D solution: {filename}")


def generate_time_series(folder, n_frames=20, nx=101, alpha=0.01):
    """Generate time series of 1D solutions for animation."""
    os.makedirs(folder, exist_ok=True)
    
    times = np.linspace(0.01, 1.0, n_frames)
    
    for i, t in enumerate(times):
        filename = os.path.join(folder, f"solution_{i:04d}.csv")
        generate_1d_solution(filename, nx=nx, t=t, alpha=alpha)
    
    print(f"Generated {n_frames} frames in {folder}")


if __name__ == '__main__':
    # Create output directory
    output_dir = "test_data"
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate test solutions
    generate_1d_solution(os.path.join(output_dir, "solution_1d.csv"))
    generate_2d_solution(os.path.join(output_dir, "solution_2d.csv"))
    generate_3d_solution(os.path.join(output_dir, "solution_3d.csv"))
    
    # Generate time series
    generate_time_series(os.path.join(output_dir, "time_series"))
    
    print("\nTest data generation complete!")
    print("\nTry visualizing with:")
    print(f"  python visualize.py {output_dir}/solution_1d.csv --dim 1")
    print(f"  python visualize.py {output_dir}/solution_2d.csv --dim 2")
    print(f"  python visualize.py {output_dir}/solution_3d.csv --dim 3 --slice 10")
