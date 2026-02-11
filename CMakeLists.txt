cmake_minimum_required(VERSION 3.14)
project(heat_equation_solver VERSION 1.0.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/heat_solver.cpp
    src/grid.cpp
    src/utils.cpp
)

# Create library
add_library(heat_solver_lib STATIC ${SOURCES})

# Main executable
add_executable(heat_solver src/main.cpp)
target_link_libraries(heat_solver heat_solver_lib)

# Examples
add_executable(example_gaussian examples/gaussian_diffusion.cpp)
target_link_libraries(example_gaussian heat_solver_lib)

# Tests
enable_testing()

add_executable(test_heat_solver tests/test_heat_solver.cpp)
target_link_libraries(test_heat_solver heat_solver_lib)
add_test(NAME HeatSolverTests COMMAND test_heat_solver)

add_executable(test_grid tests/test_grid.cpp)
target_link_libraries(test_grid heat_solver_lib)
add_test(NAME GridTests COMMAND test_grid)

# Install targets
install(TARGETS heat_solver DESTINATION bin)
install(TARGETS heat_solver_lib DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

# Print configuration
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
