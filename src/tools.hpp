#pragma once

#include <vector>
#include <chrono>
#include <random>

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"

std::vector<double> operator*(const std::vector<double>& left, double right);
std::vector<double> operator*(double left, const std::vector<double>&  right);
std::vector<double> operator+(const std::vector<double>& left, const std::vector<double>& right);
std::vector<double> operator-(const std::vector<double>& left, const std::vector<double>& right);
double operator*(const std::vector<double>& left, const std::vector<double>& right);

bool operator==(const std::vector<double>& left, const std::vector<double>& right);

double norm2(const std::vector<double>& vec);

std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec);

namespace _random {
    double getDouble(double min, double max);
    DenseMatrix getDenseMatrix(unsigned rows, unsigned cols, double density = 1, double min_el = 0, double max_el = 1, bool integer = false);
    CSRMatrix getCSRMatrix(unsigned rows, unsigned cols, double density, double min_el = 0, double max_el = 1);
    std::vector<double> getVector(unsigned rows, double min_el = 0, double max_el = 1);

    CSRMatrix getDiagonallyDominantCSRMatrix(unsigned rows, double density, double min_el = 0, double max_el = 1000);
}

namespace Jacobi {
    std::vector<double> multiply(const CSRMatrix& csr, const std::vector<double>& v);
    CSRMatrix inverseDiagonal(const CSRMatrix& csr);
}

namespace GaussSeidel {
    std::vector<double> inverseDiagonal(const CSRMatrix& csr);
}
