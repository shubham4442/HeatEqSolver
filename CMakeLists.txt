cmake_minimum_required(VERSION 3.15)
project(HeatSolver VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_STANDARD_REQUIRED ON)



add_library(SolverLib SHARED)
target_sources(SolverLib PRIVATE src/HeatSolver.cpp
                      PUBLIC  
                      FILE_SET HEADERS
                      BASE_DIRS include
                      FILES include/HeatSolver.hpp )


add_executable(HeatEqSolver)
target_sources(HeatEqSolver src/main.cpp)

target_link_libraries(HeatEqSolver PRIVATE SolverLib)