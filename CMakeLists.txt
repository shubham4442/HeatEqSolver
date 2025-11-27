cmake_minimum_required(VERSION 3.15)
project(HeatSolver VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(heat_test src/main.cpp)
