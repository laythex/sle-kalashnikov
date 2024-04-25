#pragma once

#include <vector>

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"
#include "tools.hpp"

class HessenbergMatrix {
    public:
    HessenbergMatrix(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0);
    void iterate();
    double getResidual() const;
    std::vector<double> solve() const;

    DenseMatrix getQ() const;
    DenseMatrix getR() const;

    private:
    size_t n, j;
    CSRMatrix A;
    DenseMatrix H, Q;
    std::vector<double> x0, b;
    double r0_norm;
    std::vector<std::vector<double>> V;
};