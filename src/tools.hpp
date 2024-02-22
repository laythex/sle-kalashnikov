#pragma once

#include <vector>
#include <chrono>
#include <random>

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"

double random_double(double min, double max);

DenseMatrix randomDenseMatrix(unsigned rows, unsigned cols, double density = 1, double min_el = 0, double max_el = 1, bool integer = false);
CSRMatrix randomCSRMatrix(unsigned rows, unsigned cols, double density, double min_el = 0, double max_el = 1);
std::vector<double> randomVector(unsigned rows, double min_el = 0, double max_el = 1);
