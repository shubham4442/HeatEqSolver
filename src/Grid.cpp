#include "grid.hpp"

namespace heat_equation {

//=============================================================================
// Grid1D Implementation
//=============================================================================

Grid1D::Grid1D()
    : m_dx(0.0)
    , m_x_min(0.0)
    , m_x_max(0.0)
{
}

Grid1D::Grid1D(size_t n, double x_min, double x_max)
    : m_x_min(x_min)
    , m_x_max(x_max)
{
    initialize(n, x_min, x_max);
}

Grid1D::~Grid1D()
{
}

void Grid1D::initialize(size_t n, double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
    
    // Calculate grid spacing
    m_dx = (x_max - x_min) / static_cast<double>(n - 1);
    
    // Generate coordinate array
    m_x.resize(n);
    for (size_t i = 0; i < n; ++i) {
        m_x[i] = x_min + static_cast<double>(i) * m_dx;
    }
}

size_t Grid1D::size() const
{
    return m_x.size();
}

double Grid1D::dx() const
{
    return m_dx;
}

double Grid1D::x(size_t i) const
{
    return m_x[i];
}

const std::vector<double>& Grid1D::coordinates() const
{
    return m_x;
}

//=============================================================================
// Grid2D Implementation
//=============================================================================

Grid2D::Grid2D()
    : m_dx(0.0)
    , m_dy(0.0)
{
}

Grid2D::Grid2D(size_t nx, size_t ny, double x_min, double x_max, double y_min, double y_max)
{
    initialize(nx, ny, x_min, x_max, y_min, y_max);
}

Grid2D::~Grid2D()
{
}

void Grid2D::initialize(size_t nx, size_t ny, double x_min, double x_max, double y_min, double y_max)
{
    // Calculate grid spacing
    m_dx = (x_max - x_min) / static_cast<double>(nx - 1);
    m_dy = (y_max - y_min) / static_cast<double>(ny - 1);
    
    // Generate x coordinate array
    m_x.resize(nx);
    for (size_t i = 0; i < nx; ++i) {
        m_x[i] = x_min + static_cast<double>(i) * m_dx;
    }
    
    // Generate y coordinate array
    m_y.resize(ny);
    for (size_t j = 0; j < ny; ++j) {
        m_y[j] = y_min + static_cast<double>(j) * m_dy;
    }
}

size_t Grid2D::sizeX() const
{
    return m_x.size();
}

size_t Grid2D::sizeY() const
{
    return m_y.size();
}

double Grid2D::dx() const
{
    return m_dx;
}

double Grid2D::dy() const
{
    return m_dy;
}

double Grid2D::x(size_t i) const
{
    return m_x[i];
}

double Grid2D::y(size_t j) const
{
    return m_y[j];
}

//=============================================================================
// Grid3D Implementation
//=============================================================================

Grid3D::Grid3D()
    : m_dx(0.0)
    , m_dy(0.0)
    , m_dz(0.0)
{
}

Grid3D::~Grid3D()
{
}

void Grid3D::initialize(size_t nx, size_t ny, size_t nz,
                        double x_min, double x_max,
                        double y_min, double y_max,
                        double z_min, double z_max)
{
    // Calculate grid spacing
    m_dx = (x_max - x_min) / static_cast<double>(nx - 1);
    m_dy = (y_max - y_min) / static_cast<double>(ny - 1);
    m_dz = (z_max - z_min) / static_cast<double>(nz - 1);
    
    // Generate x coordinate array
    m_x.resize(nx);
    for (size_t i = 0; i < nx; ++i) {
        m_x[i] = x_min + static_cast<double>(i) * m_dx;
    }
    
    // Generate y coordinate array
    m_y.resize(ny);
    for (size_t j = 0; j < ny; ++j) {
        m_y[j] = y_min + static_cast<double>(j) * m_dy;
    }
    
    // Generate z coordinate array
    m_z.resize(nz);
    for (size_t k = 0; k < nz; ++k) {
        m_z[k] = z_min + static_cast<double>(k) * m_dz;
    }
}

size_t Grid3D::sizeX() const
{
    return m_x.size();
}

size_t Grid3D::sizeY() const
{
    return m_y.size();
}

size_t Grid3D::sizeZ() const
{
    return m_z.size();
}

double Grid3D::dx() const
{
    return m_dx;
}

double Grid3D::dy() const
{
    return m_dy;
}

double Grid3D::dz() const
{
    return m_dz;
}

} // namespace heat_equation
