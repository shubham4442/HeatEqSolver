#!/usr/bin/env python3
"""
Heat Equation Solver Visualization Script

This script provides visualization tools for the heat equation solver output.
Supports 1D, 2D, and 3D solution visualization.

Usage:
    python visualize.py solution.csv --dim 1
    python visualize.py solution_2d.csv --dim 2
    python visualize.py solution_3d.csv --dim 3 --slice 5

Requirements:
    pip install numpy matplotlib
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.animation import FuncAnimation
import os


def load_1d_solution(filename):
    """Load 1D solution from CSV file."""
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x = data[:, 0]
    u = data[:, 1]
    return x, u


def load_2d_solution(filename):
    """Load 2D solution from CSV file."""
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x = data[:, 0]
    y = data[:, 1]
    u = data[:, 2]
    
    # Reshape to 2D grid
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    
    X = x.reshape(nx, ny)
    Y = y.reshape(nx, ny)
    U = u.reshape(nx, ny)
    
    return X, Y, U


def load_3d_solution(filename):
    """Load 3D solution from CSV file."""
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    u = data[:, 3]
    
    # Reshape to 3D grid
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    nz = len(np.unique(z))
    
    X = x.reshape(nx, ny, nz)
    Y = y.reshape(nx, ny, nz)
    Z = z.reshape(nx, ny, nz)
    U = u.reshape(nx, ny, nz)
    
    return X, Y, Z, U


def plot_1d_solution(x, u, title="1D Heat Equation Solution", save_path=None):
    """Plot 1D solution."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(x, u, 'b-', linewidth=2, label='Temperature')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('u(x)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {save_path}")
    
    plt.show()


def plot_2d_solution(X, Y, U, title="2D Heat Equation Solution", 
                     plot_type='contour', save_path=None):
    """Plot 2D solution."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Contour plot
    ax1 = axes[0]
    contour = ax1.contourf(X, Y, U, levels=50, cmap='hot')
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)
    ax1.set_title(f'{title} (Contour)', fontsize=14)
    plt.colorbar(contour, ax=ax1, label='Temperature')
    
    # 3D surface plot
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    # Remove the 2D axis we created
    axes[1].remove()
    
    surf = ax2.plot_surface(X, Y, U, cmap='hot', edgecolor='none', alpha=0.9)
    ax2.set_xlabel('x', fontsize=10)
    ax2.set_ylabel('y', fontsize=10)
    ax2.set_zlabel('u(x,y)', fontsize=10)
    ax2.set_title(f'{title} (Surface)', fontsize=14)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {save_path}")
    
    plt.show()


def plot_2d_heatmap(X, Y, U, title="2D Heat Equation Solution", save_path=None):
    """Plot 2D solution as heatmap."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(U, extent=[X.min(), X.max(), Y.min(), Y.max()],
                   origin='lower', cmap='hot', aspect='auto')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(title, fontsize=14)
    plt.colorbar(im, ax=ax, label='Temperature')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {save_path}")
    
    plt.show()


def plot_3d_slice(X, Y, Z, U, slice_index, axis='z', title="3D Heat Equation Slice",
                  save_path=None):
    """Plot a 2D slice of the 3D solution."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if axis == 'z':
        slice_data = U[:, :, slice_index]
        x_coords = X[:, :, slice_index]
        y_coords = Y[:, :, slice_index]
        z_val = Z[0, 0, slice_index]
        xlabel, ylabel = 'x', 'y'
        title_suffix = f'z = {z_val:.3f}'
    elif axis == 'y':
        slice_data = U[:, slice_index, :]
        x_coords = X[:, slice_index, :]
        y_coords = Z[:, slice_index, :]
        y_val = Y[0, slice_index, 0]
        xlabel, ylabel = 'x', 'z'
        title_suffix = f'y = {y_val:.3f}'
    else:  # axis == 'x'
        slice_data = U[slice_index, :, :]
        x_coords = Y[slice_index, :, :]
        y_coords = Z[slice_index, :, :]
        x_val = X[slice_index, 0, 0]
        xlabel, ylabel = 'y', 'z'
        title_suffix = f'x = {x_val:.3f}'
    
    im = ax.imshow(slice_data.T, origin='lower', cmap='hot', aspect='auto')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f'{title} ({title_suffix})', fontsize=14)
    plt.colorbar(im, ax=ax, label='Temperature')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {save_path}")
    
    plt.show()


def plot_3d_isosurface(X, Y, Z, U, level=0.5, title="3D Heat Equation Isosurface",
                       save_path=None):
    """Plot 3D isosurface using scatter plot approximation."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Find points close to the isosurface level
    u_min, u_max = U.min(), U.max()
    iso_value = u_min + level * (u_max - u_min)
    tolerance = 0.1 * (u_max - u_min)
    
    mask = np.abs(U - iso_value) < tolerance
    
    ax.scatter(X[mask], Y[mask], Z[mask], c=U[mask], cmap='hot', 
               alpha=0.6, s=10)
    
    ax.set_xlabel('x', fontsize=10)
    ax.set_ylabel('y', fontsize=10)
    ax.set_zlabel('z', fontsize=10)
    ax.set_title(f'{title} (level = {iso_value:.3f})', fontsize=14)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {save_path}")
    
    plt.show()


def plot_time_evolution(folder, pattern='solution_*.csv', dim=1, 
                        save_path=None):
    """Create animation of solution evolution over time."""
    import glob
    
    files = sorted(glob.glob(os.path.join(folder, pattern)))
    if not files:
        print(f"No files matching {pattern} found in {folder}")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if dim == 1:
        x, u = load_1d_solution(files[0])
        line, = ax.plot(x, u, 'b-', linewidth=2)
        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('u(x)', fontsize=12)
        ax.set_title('Heat Equation Evolution', fontsize=14)
        
        # Find global min/max for consistent y-axis
        all_u = []
        for f in files:
            _, u = load_1d_solution(f)
            all_u.extend(u)
        ax.set_ylim(min(all_u) - 0.1, max(all_u) + 0.1)
        
        def update(frame):
            x, u = load_1d_solution(files[frame])
            line.set_ydata(u)
            ax.set_title(f'Heat Equation Evolution (t = {frame})', fontsize=14)
            return line,
        
        ani = FuncAnimation(fig, update, frames=len(files), 
                           interval=100, blit=True)
        
        if save_path:
            ani.save(save_path, writer='pillow', fps=10)
            print(f"Saved animation to {save_path}")
        
        plt.show()


def compare_solutions(files, labels=None, title="Solution Comparison"):
    """Compare multiple 1D solutions on the same plot."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if labels is None:
        labels = [f"Solution {i+1}" for i in range(len(files))]
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(files)))
    
    for i, (f, label) in enumerate(zip(files, labels)):
        x, u = load_1d_solution(f)
        ax.plot(x, u, color=colors[i], linewidth=2, label=label)
    
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('u(x)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize heat equation solver output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python visualize.py solution.csv --dim 1
    python visualize.py solution_2d.csv --dim 2 --type contour
    python visualize.py solution_3d.csv --dim 3 --slice 5
    python visualize.py solution_3d.csv --dim 3 --iso 0.5
        """)
    
    parser.add_argument('filename', help='CSV file with solution data')
    parser.add_argument('--dim', type=int, default=1, choices=[1, 2, 3],
                        help='Dimension of the solution (1, 2, or 3)')
    parser.add_argument('--type', default='contour', 
                        choices=['contour', 'surface', 'heatmap'],
                        help='Plot type for 2D solutions')
    parser.add_argument('--slice', type=int, default=None,
                        help='Slice index for 3D visualization')
    parser.add_argument('--slice-axis', default='z', choices=['x', 'y', 'z'],
                        help='Axis for 3D slice')
    parser.add_argument('--iso', type=float, default=None,
                        help='Isosurface level (0-1) for 3D visualization')
    parser.add_argument('--save', type=str, default=None,
                        help='Save plot to file')
    parser.add_argument('--title', type=str, default=None,
                        help='Custom plot title')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.filename):
        print(f"Error: File {args.filename} not found")
        return
    
    title = args.title or f'{args.dim}D Heat Equation Solution'
    
    if args.dim == 1:
        x, u = load_1d_solution(args.filename)
        plot_1d_solution(x, u, title=title, save_path=args.save)
        
    elif args.dim == 2:
        X, Y, U = load_2d_solution(args.filename)
        if args.type == 'heatmap':
            plot_2d_heatmap(X, Y, U, title=title, save_path=args.save)
        else:
            plot_2d_solution(X, Y, U, title=title, save_path=args.save)
            
    elif args.dim == 3:
        X, Y, Z, U = load_3d_solution(args.filename)
        
        if args.iso is not None:
            plot_3d_isosurface(X, Y, Z, U, level=args.iso, 
                               title=title, save_path=args.save)
        elif args.slice is not None:
            plot_3d_slice(X, Y, Z, U, args.slice, axis=args.slice_axis,
                          title=title, save_path=args.save)
        else:
            # Default: show middle slice
            nz = U.shape[2]
            plot_3d_slice(X, Y, Z, U, nz // 2, axis='z',
                          title=title, save_path=args.save)


if __name__ == '__main__':
    main()
