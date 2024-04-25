#pragma once

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"
#include "HessenbergMatrix.hpp"
#include "tools.hpp"

namespace solvers {
    DenseMatrix QR(const DenseMatrix& Q, const DenseMatrix& R, const DenseMatrix& b);

    std::vector<double> fixedPointIteration(const CSRMatrix& A, const std::vector<double>& b, double tau, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> Jacobi(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> GaussSeidel(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> SymGaussSeidel(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> AccGaussSeidel(const CSRMatrix& A, const std::vector<double>& b, double rho, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> FPIAccelerated(const CSRMatrix& A, const std::vector<double>& b, double lambda_min, double lambda_max, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> GradientDescent(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> ConjugateGradient(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual);

    std::vector<double> GMRES(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual);
}