#pragma once

#include <vector>

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"
#include "tools.hpp"

// Структура для хранения поворота
struct Rotation {
    double c, s;
};

class HessenbergMatrix {
    public:
    HessenbergMatrix(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0);
    void iterate();
    double getResidual() const;
    std::vector<double> solve() const;
    
    private:
    size_t n, j;
    CSRMatrix A;
    DenseMatrix H, Q;
    std::vector<double> x0, b, q;
    double r0_norm;
    std::vector<std::vector<double>> V;
    std::vector<Rotation> rots;
};