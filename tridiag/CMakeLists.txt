cmake_minimum_required(VERSION 3.20)
project(tridiag)

find_package(GTest REQUIRED)

enable_testing()

add_executable(tridiag main.cpp tridiag_solver.hpp tridiag_solver.cpp)
target_link_libraries(tridiag GTest::gtest_main)

add_test(Alltest tridiag)
