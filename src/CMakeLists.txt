cmake_minimum_required(VERSION 3.20)
project(comparison)

add_executable(comparison comparison.cpp CSRMatrix.cpp DenseMatrix.cpp tools.cpp)
