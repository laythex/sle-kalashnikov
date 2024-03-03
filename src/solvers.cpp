#include "solvers.hpp"

DenseMatrix solvers::QR(const DenseMatrix& Q, const DenseMatrix& R, const DenseMatrix& b) {
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

std::vector<double> solvers::fixedPointIteration(const CSRMatrix& A, const std::vector<double>& b, double tau, const std::vector<double> x0, double breakpointResidual) {
    std::vector<double> x = x0;
    std::vector<double> residualVector = A * x - b;

    while (breakpointResidual < norm2(residualVector)) { 
        x = x - tau * residualVector;
        residualVector = A * x - b;  
    }

    return x;
}

std::vector<double> solvers::Jacobi(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double> x0, double breakpointResidual) {
    std::vector<double> x = x0;
    CSRMatrix d = JacobiTools::inverseDiagonal(A);
    double residual = norm2(A * x - b);

    while (breakpointResidual < residual) {
        x = d * (b - JacobiTools::multiply(A, x));
        residual = norm2(A * x - b);
    }

    return x;
}

std::vector<double> solvers::GaussSeidel(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double> x0, double breakpointResidual) {
    std::vector<double> x = x0;
    std::vector<double> d = GaussSeidelTools::inverseDiagonal(A);
    double residual = norm2(A * x - b), tmp;

    while (breakpointResidual < residual) {
        for (unsigned i = 0; i < x.size(); i++) {
            tmp = 0;
            for (unsigned j = A.rowsAt(i); j < A.rowsAt(i + 1); j++) {
                if (i != A.colsAt(j)) {
                    tmp += A.valsAt(j) * x[A.colsAt(j)];
                }
            }
            x[i] = d[i] * (b[i] - tmp);
        }
        residual = norm2(A * x - b);
    }

    return x;
}

std::vector<double> solvers::FPIAccelerated(const CSRMatrix& A, const std::vector<double>& b, double lambda, const std::vector<double>& roots, const std::vector<size_t>& permutations, const std::vector<double> x0, double breakpointResidual) {
    std::vector<double> x = x0;
    std::vector<double> residualVector = A * x - b;
    size_t n = roots.size();
    
    double tau;
    while (breakpointResidual < norm2(residualVector)) { 
        for (size_t i = 0; i < n; i++) {
            tau = 1 / (lambda / 2 * (1 + roots[permutations[i]]));
            x = x - tau * residualVector;
            residualVector = A * x - b;
        }
    }

    return x;
}

