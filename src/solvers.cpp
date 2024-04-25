#include "solvers.hpp"

DenseMatrix solvers::QR(const DenseMatrix& Q, const DenseMatrix& R, const DenseMatrix& b) {
    DenseMatrix v = Q.transpose() * b;
    unsigned n = b.getRows();
    std::vector<double> x(n);
    double tmp;

    for (unsigned i = n; i --> 0;) {
        tmp = 0;
        for (unsigned j = n - 1; j > i; j--) {
            tmp += R(i, j) * x[j];
        }
        x[i] = (v(i, 0) - tmp) / R(i, i); 
    }

    return DenseMatrix(x, 1);
}

std::vector<double> solvers::fixedPointIteration(const CSRMatrix& A, const std::vector<double>& b, double tau, const std::vector<double>& x0, double breakpointResidual) {
    std::vector<double> x = x0;
    std::vector<double> residualVector = A * x - b;

    while (breakpointResidual < norm2(residualVector)) { 
        x = x - tau * residualVector;
        residualVector = A * x - b;  
    }

    return x;
}

std::vector<double> solvers::Jacobi(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual) {
    std::vector<double> x = x0;
    CSRMatrix d = JacobiTools::inverseDiagonal(A);
    double residual = norm2(A * x - b);

    while (breakpointResidual < residual) {
        x = d * (b - JacobiTools::multiply(A, x));
        residual = norm2(A * x - b);
    }

    return x;
}

std::vector<double> solvers::GaussSeidel(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual) {
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

std::vector<double> solvers::SymGaussSeidel(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual) {
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
 
        for (unsigned i = x.size(); i --> 0;) {
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
std::vector<double> solvers::AccGaussSeidel(const CSRMatrix& A, const std::vector<double>& b, double rho, const std::vector<double>& x0, double breakpointResidual) {
    std::vector<double> xp = x0, x, xn;
    std::vector<double> d = GaussSeidelTools::inverseDiagonal(A);
    double residual, mp = 1, m = 1 / rho, mn;

    x = AccGaussSeidelTools::iterate(A, b, d, xp);
    residual = norm2(A * x - b);

    while (breakpointResidual < residual) {
        mn = 2 / rho * m - mp;
        xn = AccGaussSeidelTools::iterate(A, b, d, x);
        xn = 2 / rho * m / mn * xn - mp / mn * xp;

        xp = x;
        x = xn;
        mp = m;
        m = mn;

        residual = norm2(A * xn - b);
    }

    return xn;
}

std::vector<double> solvers::FPIAccelerated(const CSRMatrix& A, const std::vector<double>& b, double lambda_min, double lambda_max, const std::vector<double>& x0, double breakpointResidual) {
    std::vector<double> x = x0;
    std::vector<double> residualVector = A * x - b;
    
    double tau;
    double k1 = (lambda_max + lambda_min) / 2, k2 = (lambda_max - lambda_min) / 2;
    
    while (breakpointResidual < norm2(residualVector)) { 
        for (size_t i = 0; i < FPIATools::nroots; i++) {
            tau = 1 / (k1 + k2 * FPIATools::roots[FPIATools::permutations[i]]);
            x = x - tau * residualVector;
            residualVector = A * x - b;
        }
    }

    return x;
}

std::vector<double> solvers::GradientDescent(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual) {
    std::vector<double> x = x0;
    std::vector<double> residualVector = A * x - b;
    double tau;

    while (breakpointResidual < norm2(residualVector)) { 
        tau = residualVector * residualVector / (residualVector * (A * residualVector));
        x = x - tau * residualVector;
        residualVector = A * x - b;
    }

    return x;
}

std::vector<double> solvers::ConjugateGradient(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual) {
    std::vector<double> x = x0;
    std::vector<double> r = A * x - b, r0;
    std::vector<double> d = r;
    double alpha, beta;

    while (breakpointResidual < norm2(r)) { 
        alpha = r * r / (d * (A * d));
        x = x - alpha * d;
        r0 = r;
        r = A * x - b;
        beta = r * r / (r0 * r0);
        d = r + beta * d;
    }

    return x;
}

std::vector<double> solvers::GMRES(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0, double breakpointResidual) {
    HessenbergMatrix hm(A, b, x0);

    while (breakpointResidual < hm.getResidual()) {
        hm.iterate();
    }

    return hm.solve();
}

