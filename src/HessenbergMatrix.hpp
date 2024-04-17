#pragma once

#include <vector>

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"
#include "tools.hpp"

class HessenbergMatrix {
    public:
    HessenbergMatrix(const CSRMatrix& A, const std::vector<double>& r0);
    void iterate();

    DenseMatrix H;
    
    private:
    size_t n, j;
    CSRMatrix A;
    std::vector<double> r0;
    std::vector<std::vector<double>> V;
    std::vector<DenseMatrix> rots;
};