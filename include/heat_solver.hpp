#ifndef HEAT_SOLVER_HPP
#define HEAT_SOLVER_HPP

#include <vector>
#include <string>
#include "grid.hpp"

namespace heat_equation {

/**
 * @brief Boundary condition types for the heat equation
 */
enum class BoundaryCondition {
    Dirichlet,
    Neumann,
    Periodic
};

/**
 * @brief Configuration parameters for the heat equation solver
 */
struct SolverConfig {
    double alpha;           // Thermal diffusivity
    double dx;              // Spatial step size
    double dy;              // Spatial step size (for 2D/3D)
    double dz;              // Spatial step size (for 3D)
    double dt;              // Time step size
    size_t nx;              // Number of spatial grid points
    size_t ny;              // Number of spatial grid points (for 2D)
    size_t nz;              // Number of spatial grid points (for 3D)
    double t_final;         // Final simulation time
    BoundaryCondition bc;   // Boundary condition type
    
    // Domain bounds (optional, used for grid generation)
    double x_min = 0.0;
    double x_max = 1.0;
    double y_min = 0.0;
    double y_max = 1.0;
    double z_min = 0.0;
    double z_max = 1.0;
};

/**
 * @brief Abstract base class for heat equation solvers
 */
class HeatSolverBase {
public:
    HeatSolverBase();
    virtual ~HeatSolverBase();

    virtual void initialize(const SolverConfig& config) = 0;
    virtual void setInitialCondition(const std::vector<double>& u0) = 0;
    virtual void setBoundaryConditions(double left, double right) = 0;
    virtual void solve() = 0;
    virtual std::vector<double> getSolution() const = 0;
    virtual void exportSolution(const std::string& filename) const = 0;
    
    double getCurrentTime() const { return m_time; }
    const SolverConfig& getConfig() const { return m_config; }

protected:
    SolverConfig m_config;
    bool m_initialized;
    double m_time;
};

/**
 * @brief 1D Heat equation solver using explicit finite difference method
 */
class HeatSolver1D : public HeatSolverBase {
public:
    HeatSolver1D();
    ~HeatSolver1D() override;

    void initialize(const SolverConfig& config) override;
    void setInitialCondition(const std::vector<double>& u0) override;
    void setBoundaryConditions(double left, double right) override;
    void solve() override;
    std::vector<double> getSolution() const override;
    void exportSolution(const std::string& filename) const override;

    void step();
    double computeStabilityCondition() const;
    
    const Grid1D& getGrid() const { return m_grid; }

private:
    Grid1D m_grid;
    std::vector<double> m_u;
    std::vector<double> m_u_new;
    double m_bc_left;
    double m_bc_right;
};

/**
 * @brief 2D Heat equation solver using explicit finite difference method
 */
class HeatSolver2D : public HeatSolverBase {
public:
    HeatSolver2D();
    ~HeatSolver2D() override;

    void initialize(const SolverConfig& config) override;
    void setInitialCondition(const std::vector<double>& u0) override;
    void setBoundaryConditions(double left, double right) override;
    void solve() override;
    std::vector<double> getSolution() const override;
    void exportSolution(const std::string& filename) const override;

    void step();
    double computeStabilityCondition() const;
    
    const Grid2D& getGrid() const { return m_grid; }

private:
    Grid2D m_grid;
    std::vector<std::vector<double>> m_u;
    std::vector<std::vector<double>> m_u_new;
};

/**
 * @brief 3D Heat equation solver using explicit finite difference method
 */
class HeatSolver3D : public HeatSolverBase {
public:
    HeatSolver3D();
    ~HeatSolver3D() override;

    void initialize(const SolverConfig& config) override;
    void setInitialCondition(const std::vector<double>& u0) override;
    void setBoundaryConditions(double left, double right) override;
    void solve() override;
    std::vector<double> getSolution() const override;
    void exportSolution(const std::string& filename) const override;

    void step();
    double computeStabilityCondition() const;
    
    const Grid3D& getGrid() const { return m_grid; }
    
    // Export a 2D slice at given z-index
    void exportSlice(const std::string& filename, size_t z_index) const;

private:
    Grid3D m_grid;
    std::vector<std::vector<std::vector<double>>> m_u;
    std::vector<std::vector<std::vector<double>>> m_u_new;
};

} // namespace heat_equation

#endif // HEAT_SOLVER_HPP
