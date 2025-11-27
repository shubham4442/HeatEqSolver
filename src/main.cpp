#include <iostream>
#include <vector>
#include <cmath>

int main() {
    std::cout << "C++ environment OK! \n";

    // quick runtime test
    std::vector<double> x(5);
    for (int i = 0; i < 5; i++) x[i] = std::sin(i);

    std::cout << "Vector: ";
    for (double v : x) std::cout << v << " ";
    std::cout << "\n";

    return 0;
}
