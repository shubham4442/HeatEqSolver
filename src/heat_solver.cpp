#include "heat_solver.hpp"
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace heat_equation {

//=============================================================================
// HeatSolverBase Implementation
//=============================================================================

HeatSolverBase::HeatSolverBase()
    : m_initialized(false)
    , m_time(0.0)
{
}

HeatSolverBase::~HeatSolverBase()
{
}

//=============================================================================
// HeatSolver1D Implementation
//=============================================================================

HeatSolver1D::HeatSolver1D()
    : HeatSolverBase()
    , m_bc_left(0.0)
    , m_bc_right(0.0)
{
}

HeatSolver1D::~HeatSolver1D()
{
}

void HeatSolver1D::initialize(const SolverConfig& config)
{
    m_config = config;
    
    // Initialize the grid
    m_grid.initialize(config.nx, config.x_min, config.x_max);
    
    // Allocate memory for solution vectors
    m_u.resize(config.nx, 0.0);
    m_u_new.resize(config.nx, 0.0);
    
    // Check stability condition
    double r = computeStabilityCondition();
    if (r > 0.5) {
        std::cerr << "Warning: Stability condition violated! r = " << r 
                  << " (should be <= 0.5)" << std::endl;
    }
    
    m_time = 0.0;
    m_initialized = true;
}

void HeatSolver1D::setInitialCondition(const std::vector<double>& u0)
{
    if (u0.size() != m_config.nx) {
        throw std::invalid_argument("Initial condition size does not match grid size");
    }
    m_u = u0;
    m_u_new = u0;
}

void HeatSolver1D::setBoundaryConditions(double left, double right)
{
    m_bc_left = left;
    m_bc_right = right;
}

void HeatSolver1D::solve()
{
    if (!m_initialized) {
        throw std::runtime_error("Solver not initialized");
    }
    
    size_t num_steps = static_cast<size_t>(m_config.t_final / m_config.dt);
    
    for (size_t n = 0; n < num_steps; ++n) {
        step();
        m_time += m_config.dt;
    }
}

void HeatSolver1D::step()
{
    // Stability parameter r = alpha * dt / dx^2
    double r = m_config.alpha * m_config.dt / (m_config.dx * m_config.dx);
    
    // Apply explicit finite difference scheme for interior points
    for (size_t i = 1; i < m_config.nx - 1; ++i) {
        m_u_new[i] = m_u[i] + r * (m_u[i + 1] - 2.0 * m_u[i] + m_u[i - 1]);
    }
    
    // Apply boundary conditions
    if (m_config.bc == BoundaryCondition::Dirichlet) {
        m_u_new[0] = m_bc_left;
        m_u_new[m_config.nx - 1] = m_bc_right;
    } else if (m_config.bc == BoundaryCondition::Neumann) {
        m_u_new[0] = m_u_new[1] - m_bc_left * m_config.dx;
        m_u_new[m_config.nx - 1] = m_u_new[m_config.nx - 2] + m_bc_right * m_config.dx;
    } else if (m_config.bc == BoundaryCondition::Periodic) {
        m_u_new[0] = m_u[0] + r * (m_u[1] - 2.0 * m_u[0] + m_u[m_config.nx - 2]);
        m_u_new[m_config.nx - 1] = m_u_new[0];
    }
    
    std::swap(m_u, m_u_new);
}

std::vector<double> HeatSolver1D::getSolution() const
{
    return m_u;
}

void HeatSolver1D::exportSolution(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    file << "x,u" << std::endl;
    for (size_t i = 0; i < m_config.nx; ++i) {
        file << m_grid.x(i) << "," << m_u[i] << std::endl;
    }
    
    file.close();
}

double HeatSolver1D::computeStabilityCondition() const
{
    return m_config.alpha * m_config.dt / (m_config.dx * m_config.dx);
}

//=============================================================================
// HeatSolver2D Implementation
//=============================================================================

HeatSolver2D::HeatSolver2D()
    : HeatSolverBase()
{
}

HeatSolver2D::~HeatSolver2D()
{
}

void HeatSolver2D::initialize(const SolverConfig& config)
{
    m_config = config;
    
    // Initialize the grid
    m_grid.initialize(config.nx, config.ny, 
                      config.x_min, config.x_max,
                      config.y_min, config.y_max);
    
    // Allocate 2D arrays
    m_u.resize(config.nx, std::vector<double>(config.ny, 0.0));
    m_u_new.resize(config.nx, std::vector<double>(config.ny, 0.0));
    
    // Check stability condition for 2D (r <= 0.25)
    double r = computeStabilityCondition();
    if (r > 0.25) {
        std::cerr << "Warning: 2D Stability condition violated! r = " << r 
                  << " (should be <= 0.25)" << std::endl;
    }
    
    m_time = 0.0;
    m_initialized = true;
}

void HeatSolver2D::setInitialCondition(const std::vector<double>& u0)
{
    if (u0.size() != m_config.nx * m_config.ny) {
        throw std::invalid_argument("Initial condition size does not match grid size");
    }
    
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            m_u[i][j] = u0[i * m_config.ny + j];
            m_u_new[i][j] = u0[i * m_config.ny + j];
        }
    }
}

void HeatSolver2D::setBoundaryConditions(double left, double right)
{
    for (size_t i = 0; i < m_config.nx; ++i) {
        m_u[i][0] = left;
        m_u[i][m_config.ny - 1] = right;
        m_u_new[i][0] = left;
        m_u_new[i][m_config.ny - 1] = right;
    }
    for (size_t j = 0; j < m_config.ny; ++j) {
        m_u[0][j] = left;
        m_u[m_config.nx - 1][j] = right;
        m_u_new[0][j] = left;
        m_u_new[m_config.nx - 1][j] = right;
    }
}

void HeatSolver2D::solve()
{
    if (!m_initialized) {
        throw std::runtime_error("Solver not initialized");
    }
    
    size_t num_steps = static_cast<size_t>(m_config.t_final / m_config.dt);
    
    for (size_t n = 0; n < num_steps; ++n) {
        step();
        m_time += m_config.dt;
    }
}

void HeatSolver2D::step()
{
    double dx = m_grid.dx();
    double dy = m_grid.dy();
    double rx = m_config.alpha * m_config.dt / (dx * dx);
    double ry = m_config.alpha * m_config.dt / (dy * dy);
    
    for (size_t i = 1; i < m_config.nx - 1; ++i) {
        for (size_t j = 1; j < m_config.ny - 1; ++j) {
            m_u_new[i][j] = m_u[i][j] 
                + rx * (m_u[i + 1][j] - 2.0 * m_u[i][j] + m_u[i - 1][j])
                + ry * (m_u[i][j + 1] - 2.0 * m_u[i][j] + m_u[i][j - 1]);
        }
    }
    
    std::swap(m_u, m_u_new);
}

std::vector<double> HeatSolver2D::getSolution() const
{
    std::vector<double> result(m_config.nx * m_config.ny);
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            result[i * m_config.ny + j] = m_u[i][j];
        }
    }
    return result;
}

void HeatSolver2D::exportSolution(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // Export with coordinates
    file << "x,y,u" << std::endl;
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            file << m_grid.x(i) << "," << m_grid.y(j) << "," << m_u[i][j] << std::endl;
        }
    }
    
    file.close();
}

double HeatSolver2D::computeStabilityCondition() const
{
    double dx = m_config.dx > 0 ? m_config.dx : (m_config.x_max - m_config.x_min) / (m_config.nx - 1);
    double dy = m_config.dy > 0 ? m_config.dy : (m_config.y_max - m_config.y_min) / (m_config.ny - 1);
    return m_config.alpha * m_config.dt * (1.0/(dx*dx) + 1.0/(dy*dy));
}

//=============================================================================
// HeatSolver3D Implementation
//=============================================================================

HeatSolver3D::HeatSolver3D()
    : HeatSolverBase()
{
}

HeatSolver3D::~HeatSolver3D()
{
}

void HeatSolver3D::initialize(const SolverConfig& config)
{
    m_config = config;
    
    // Initialize the 3D grid
    m_grid.initialize(config.nx, config.ny, config.nz,
                      config.x_min, config.x_max,
                      config.y_min, config.y_max,
                      config.z_min, config.z_max);
    
    // Allocate 3D arrays
    m_u.resize(config.nx, 
        std::vector<std::vector<double>>(config.ny, 
            std::vector<double>(config.nz, 0.0)));
    m_u_new.resize(config.nx, 
        std::vector<std::vector<double>>(config.ny, 
            std::vector<double>(config.nz, 0.0)));
    
    // Check stability condition for 3D (r <= 1/6)
    double r = computeStabilityCondition();
    if (r > 1.0/6.0) {
        std::cerr << "Warning: 3D Stability condition violated! r = " << r 
                  << " (should be <= " << 1.0/6.0 << ")" << std::endl;
    }
    
    m_time = 0.0;
    m_initialized = true;
}

void HeatSolver3D::setInitialCondition(const std::vector<double>& u0)
{
    if (u0.size() != m_config.nx * m_config.ny * m_config.nz) {
        throw std::invalid_argument("Initial condition size does not match grid size");
    }
    
    // Reshape 1D input to 3D grid
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            for (size_t k = 0; k < m_config.nz; ++k) {
                size_t idx = i * m_config.ny * m_config.nz + j * m_config.nz + k;
                m_u[i][j][k] = u0[idx];
                m_u_new[i][j][k] = u0[idx];
            }
        }
    }
}

void HeatSolver3D::setBoundaryConditions(double left, double right)
{
    // Set all 6 faces of the 3D domain
    // X boundaries
    for (size_t j = 0; j < m_config.ny; ++j) {
        for (size_t k = 0; k < m_config.nz; ++k) {
            m_u[0][j][k] = left;
            m_u[m_config.nx - 1][j][k] = right;
            m_u_new[0][j][k] = left;
            m_u_new[m_config.nx - 1][j][k] = right;
        }
    }
    // Y boundaries
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t k = 0; k < m_config.nz; ++k) {
            m_u[i][0][k] = left;
            m_u[i][m_config.ny - 1][k] = right;
            m_u_new[i][0][k] = left;
            m_u_new[i][m_config.ny - 1][k] = right;
        }
    }
    // Z boundaries
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            m_u[i][j][0] = left;
            m_u[i][j][m_config.nz - 1] = right;
            m_u_new[i][j][0] = left;
            m_u_new[i][j][m_config.nz - 1] = right;
        }
    }
}

void HeatSolver3D::solve()
{
    if (!m_initialized) {
        throw std::runtime_error("Solver not initialized");
    }
    
    size_t num_steps = static_cast<size_t>(m_config.t_final / m_config.dt);
    
    for (size_t n = 0; n < num_steps; ++n) {
        step();
        m_time += m_config.dt;
    }
}

void HeatSolver3D::step()
{
    double dx = m_grid.dx();
    double dy = m_grid.dy();
    double dz = m_grid.dz();
    
    double rx = m_config.alpha * m_config.dt / (dx * dx);
    double ry = m_config.alpha * m_config.dt / (dy * dy);
    double rz = m_config.alpha * m_config.dt / (dz * dz);
    
    // Apply 3D explicit finite difference scheme for interior points
    for (size_t i = 1; i < m_config.nx - 1; ++i) {
        for (size_t j = 1; j < m_config.ny - 1; ++j) {
            for (size_t k = 1; k < m_config.nz - 1; ++k) {
                m_u_new[i][j][k] = m_u[i][j][k]
                    + rx * (m_u[i+1][j][k] - 2.0*m_u[i][j][k] + m_u[i-1][j][k])
                    + ry * (m_u[i][j+1][k] - 2.0*m_u[i][j][k] + m_u[i][j-1][k])
                    + rz * (m_u[i][j][k+1] - 2.0*m_u[i][j][k] + m_u[i][j][k-1]);
            }
        }
    }
    
    std::swap(m_u, m_u_new);
}

std::vector<double> HeatSolver3D::getSolution() const
{
    std::vector<double> result(m_config.nx * m_config.ny * m_config.nz);
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            for (size_t k = 0; k < m_config.nz; ++k) {
                size_t idx = i * m_config.ny * m_config.nz + j * m_config.nz + k;
                result[idx] = m_u[i][j][k];
            }
        }
    }
    return result;
}

void HeatSolver3D::exportSolution(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // Export with coordinates (can be large!)
    file << "x,y,z,u" << std::endl;
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            for (size_t k = 0; k < m_config.nz; ++k) {
                file << m_grid.dx() * i << "," 
                     << m_grid.dy() * j << "," 
                     << m_grid.dz() * k << ","
                     << m_u[i][j][k] << std::endl;
            }
        }
    }
    
    file.close();
}

void HeatSolver3D::exportSlice(const std::string& filename, size_t z_index) const
{
    if (z_index >= m_config.nz) {
        throw std::out_of_range("z_index out of range");
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // Export a 2D slice at given z-index
    file << "x,y,u" << std::endl;
    for (size_t i = 0; i < m_config.nx; ++i) {
        for (size_t j = 0; j < m_config.ny; ++j) {
            file << m_grid.dx() * i << "," 
                 << m_grid.dy() * j << ","
                 << m_u[i][j][z_index] << std::endl;
        }
    }
    
    file.close();
}

double HeatSolver3D::computeStabilityCondition() const
{
    double dx = m_config.dx > 0 ? m_config.dx : (m_config.x_max - m_config.x_min) / (m_config.nx - 1);
    double dy = m_config.dy > 0 ? m_config.dy : (m_config.y_max - m_config.y_min) / (m_config.ny - 1);
    double dz = m_config.dz > 0 ? m_config.dz : (m_config.z_max - m_config.z_min) / (m_config.nz - 1);
    return m_config.alpha * m_config.dt * (1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz));
}

} // namespace heat_equation
