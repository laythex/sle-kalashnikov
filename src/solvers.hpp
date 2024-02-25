#pragma once

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"
#include "tools.hpp"

namespace solvers {

    DenseMatrix QR(const DenseMatrix& Q, const DenseMatrix& R, const DenseMatrix& b);

    std::vector<double> fixedPointIteration(const CSRMatrix& A, const std::vector<double>& b, double tau, const std::vector<double> x0, double breakpointResidual);

    std::vector<double> Jacobi(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double> x0, double breakpointResidual);

    std::vector<double> GaussSeidel(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double> x0, double breakpointResidual);

}