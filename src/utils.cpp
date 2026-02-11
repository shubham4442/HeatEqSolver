#include "utils.hpp"
#include <fstream>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace heat_equation {
namespace utils {

void writeToCsv(const std::string& filename,
                const std::vector<double>& x,
                const std::vector<double>& u)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    file << "x,u" << std::endl;
    size_t n = std::min(x.size(), u.size());
    for (size_t i = 0; i < n; ++i) {
        file << x[i] << "," << u[i] << std::endl;
    }
    
    file.close();
}

void writeToCsv2D(const std::string& filename,
                  const std::vector<std::vector<double>>& u)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    for (size_t i = 0; i < u.size(); ++i) {
        for (size_t j = 0; j < u[i].size(); ++j) {
            file << u[i][j];
            if (j < u[i].size() - 1) file << ",";
        }
        file << std::endl;
    }
    
    file.close();
}

std::vector<double> readFromFile(const std::string& filename)
{
    std::vector<double> values;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double value;
        while (iss >> value) {
            values.push_back(value);
            if (iss.peek() == ',' || iss.peek() == ' ') {
                iss.ignore();
            }
        }
    }
    
    file.close();
    return values;
}

std::vector<double> generateInitialCondition(
    size_t n,
    double x_min,
    double x_max,
    std::function<double(double)> func)
{
    std::vector<double> u0(n);
    double dx = (x_max - x_min) / static_cast<double>(n - 1);
    
    for (size_t i = 0; i < n; ++i) {
        double x = x_min + static_cast<double>(i) * dx;
        u0[i] = func(x);
    }
    
    return u0;
}

double computeL2Error(const std::vector<double>& u,
                      const std::vector<double>& u_exact)
{
    if (u.size() != u_exact.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < u.size(); ++i) {
        double diff = u[i] - u_exact[i];
        sum += diff * diff;
    }
    
    return std::sqrt(sum / static_cast<double>(u.size()));
}

double computeMaxError(const std::vector<double>& u,
                       const std::vector<double>& u_exact)
{
    if (u.size() != u_exact.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    
    double max_error = 0.0;
    for (size_t i = 0; i < u.size(); ++i) {
        double error = std::abs(u[i] - u_exact[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    
    return max_error;
}

//=============================================================================
// Timer Implementation
//=============================================================================

Timer::Timer()
    : m_start(0.0)
    , m_elapsed(0.0)
    , m_running(false)
{
}

Timer::~Timer()
{
}

void Timer::start()
{
    auto now = std::chrono::high_resolution_clock::now();
    m_start = std::chrono::duration<double>(now.time_since_epoch()).count();
    m_running = true;
}

void Timer::stop()
{
    if (m_running) {
        auto now = std::chrono::high_resolution_clock::now();
        double end = std::chrono::duration<double>(now.time_since_epoch()).count();
        m_elapsed += end - m_start;
        m_running = false;
    }
}

double Timer::elapsed() const
{
    if (m_running) {
        auto now = std::chrono::high_resolution_clock::now();
        double current = std::chrono::duration<double>(now.time_since_epoch()).count();
        return m_elapsed + (current - m_start);
    }
    return m_elapsed;
}

void Timer::reset()
{
    m_elapsed = 0.0;
    m_running = false;
    m_start = 0.0;
}

} // namespace utils
} // namespace heat_equation
