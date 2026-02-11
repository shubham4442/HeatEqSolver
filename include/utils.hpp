#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include <functional>

namespace heat_equation {
namespace utils {

/**
 * @brief Write solution data to CSV file
 */
void writeToCsv(const std::string& filename,
                const std::vector<double>& x,
                const std::vector<double>& u);

/**
 * @brief Write 2D solution data to CSV file
 */
void writeToCsv2D(const std::string& filename,
                  const std::vector<std::vector<double>>& u);

/**
 * @brief Read initial condition from file
 */
std::vector<double> readFromFile(const std::string& filename);

/**
 * @brief Generate initial condition using a function
 */
std::vector<double> generateInitialCondition(
    size_t n,
    double x_min,
    double x_max,
    std::function<double(double)> func);

/**
 * @brief Compute L2 norm of the difference between two vectors
 */
double computeL2Error(const std::vector<double>& u,
                      const std::vector<double>& u_exact);

/**
 * @brief Compute maximum norm of the difference between two vectors
 */
double computeMaxError(const std::vector<double>& u,
                       const std::vector<double>& u_exact);

/**
 * @brief Timer class for performance measurement
 */
class Timer {
public:
    Timer();
    ~Timer();

    void start();
    void stop();
    double elapsed() const;
    void reset();

private:
    double m_start;
    double m_elapsed;
    bool m_running;
};

} // namespace utils
} // namespace heat_equation

#endif // UTILS_HPP
