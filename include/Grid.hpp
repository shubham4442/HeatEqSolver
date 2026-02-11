#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <cstddef>

namespace heat_equation {

/**
 * @brief 1D computational grid
 */
class Grid1D {
public:
    Grid1D();
    Grid1D(size_t n, double x_min, double x_max);
    ~Grid1D();

    void initialize(size_t n, double x_min, double x_max);
    
    size_t size() const;
    double dx() const;
    double x(size_t i) const;
    
    const std::vector<double>& coordinates() const;

private:
    std::vector<double> m_x;
    double m_dx;
    double m_x_min;
    double m_x_max;
};

/**
 * @brief 2D computational grid
 */
class Grid2D {
public:
    Grid2D();
    Grid2D(size_t nx, size_t ny, double x_min, double x_max, double y_min, double y_max);
    ~Grid2D();

    void initialize(size_t nx, size_t ny, double x_min, double x_max, double y_min, double y_max);
    
    size_t sizeX() const;
    size_t sizeY() const;
    double dx() const;
    double dy() const;
    double x(size_t i) const;
    double y(size_t j) const;

private:
    std::vector<double> m_x;
    std::vector<double> m_y;
    double m_dx;
    double m_dy;
};

/**
 * @brief 3D computational grid
 */
class Grid3D {
public:
    Grid3D();
    ~Grid3D();

    void initialize(size_t nx, size_t ny, size_t nz,
                    double x_min, double x_max,
                    double y_min, double y_max,
                    double z_min, double z_max);

    size_t sizeX() const;
    size_t sizeY() const;
    size_t sizeZ() const;
    double dx() const;
    double dy() const;
    double dz() const;

private:
    std::vector<double> m_x;
    std::vector<double> m_y;
    std::vector<double> m_z;
    double m_dx;
    double m_dy;
    double m_dz;
};

} // namespace heat_equation

#endif // GRID_HPP
