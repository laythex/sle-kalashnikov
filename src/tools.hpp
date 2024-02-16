#pragma once

#include <vector>
#include <chrono>
#include <random>

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"

double random_double(double min, double max);

DenseMatrix generateRandomDenseMatrix(unsigned rows, unsigned cols, double density = 1, double min_el = 0, double max_el = 1);
CSRMatrix generateRandomCSRMatrix(unsigned rows, unsigned cols, double density, double min_el = 0, double max_el = 1);
std::vector<double> generateRandomVector(unsigned rows, double min_el = 0, double max_el = 1);