#include "solvers.hpp"

DenseMatrix solveQR(const DenseMatrix& Q, const DenseMatrix& R, const DenseMatrix& b) {
    DenseMatrix v = Q.transpose() * b;
    unsigned n = b.getRows();
    std::vector<double> x(n);
    double tmp;
    for (int i = n - 1; i >= 0; i--) {
        tmp = 0;
        for (unsigned j = n - 1; j > i; j--) {
            tmp += R(i, j) * x[j];
        }
        x[i] = (v(i, 0) - tmp) / R(i, i); 
    }
    return DenseMatrix(x, 1);
}