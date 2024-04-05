#include <gtest/gtest.h>

#include "CSRMatrix.hpp"
#include "tools.hpp"
#include "solvers.hpp"
#include "HouseholderQR.hpp"

TEST(solve, QR) {
    DenseMatrix A = DenseMatrix({1, 3, 4,
                                 5, 3, 2,
                                 2, 6, 1}, 3);
    DenseMatrix x_real = DenseMatrix({1, 1, 1}, 1);
    DenseMatrix b = A * x_real;

    HouseholderQR hhqr = HouseholderQR(A);
    DenseMatrix x_calc = solvers::QR(hhqr.getQ(), hhqr.getR(), b);
	
    EXPECT_TRUE(x_calc == x_real);
}

TEST(solve, FPI) {
    unsigned n = _random::getUnsigned(1, 1e2);
    CSRMatrix A = _random::getTestMatrix(n);
	std::vector<double> x_real = _random::getVector(n);
    std::vector<double> b = A * x_real;
    std::vector<double> x0 = _random::getVector(n);

    std::vector<double> x_calc = solvers::fixedPointIteration(A, b, 0.1, x0, 1e-14);

    EXPECT_TRUE(x_calc == x_real);
}

TEST(solve, Jacobi) {
    unsigned n = _random::getUnsigned(1, 1e2);
    CSRMatrix A = _random::getTestMatrix(n);
	std::vector<double> x_real = _random::getVector(n);
    std::vector<double> b = A * x_real;
    std::vector<double> x0 = _random::getVector(n);

    std::vector<double> x_calc = solvers::Jacobi(A, b, x0, 1e-14);

    EXPECT_TRUE(x_calc == x_real);
}

TEST(solve, GaussSeidel) {
    unsigned n = _random::getUnsigned(1, 1e2);
    CSRMatrix A = _random::getTestMatrix(n);
	std::vector<double> x_real = _random::getVector(n);
    std::vector<double> b = A * x_real;
    std::vector<double> x0 = _random::getVector(n);

    std::vector<double> x_calc = solvers::GaussSeidel(A, b, x0, 1e-14);

    EXPECT_TRUE(x_calc == x_real);
}

TEST(solve, SymGaussSeidel) {
    unsigned n = _random::getUnsigned(1, 1e2);
    CSRMatrix A = _random::getTestMatrix(n);
	std::vector<double> x_real = _random::getVector(n);
    std::vector<double> b = A * x_real;
    std::vector<double> x0 = _random::getVector(n);

    std::vector<double> x_calc = solvers::SymGaussSeidel(A, b, x0, 1e-14);

    EXPECT_TRUE(x_calc == x_real);
}

TEST(solve, AccGaussSeidel) {
    unsigned n = _random::getUnsigned(1, 1e3);
    CSRMatrix A = _random::getTestMatrix(n);
	std::vector<double> x_real = _random::getVector(n);
    std::vector<double> b = A * x_real;
    std::vector<double> x0 = _random::getVector(n);

    double rho = 0.9;

    std::vector<double> x_calc = solvers::AccGaussSeidel(A, b, rho, x0, 1e-14);
    EXPECT_TRUE(x_calc == x_real);
}

TEST(solve, FPIA) {
    unsigned n = _random::getUnsigned(1, 1e2);
    CSRMatrix A = _random::getTestMatrix(n);
	std::vector<double> x_real = _random::getVector(n);
    std::vector<double> b = A * x_real;
    std::vector<double> x0 = _random::getVector(n);

    double lambda = calcMaxEigenvalue(A, 1e-5);

    std::vector<double> x_calc = solvers::FPIAccelerated(A, b, 0, lambda, x0, 1e-14);

    EXPECT_TRUE(x_calc == x_real);
}

TEST(solve, GradientDescent) {
    unsigned n = _random::getUnsigned(1, 1e2);
    CSRMatrix A = _random::getTestMatrix(n);
	std::vector<double> x_real = _random::getVector(n);
    std::vector<double> b = A * x_real;
    std::vector<double> x0 = _random::getVector(n);

    std::vector<double> x_calc = solvers::GradientDescent(A, b, x0, 1e-14);

    EXPECT_TRUE(x_calc == x_real);
}