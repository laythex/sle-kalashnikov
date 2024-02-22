#include <gtest/gtest.h>

#include "DenseMatrix.hpp"
#include "HouseholderQR.hpp"
#include "solvers.hpp"

TEST(QRsolver, main) {
	DenseMatrix A = DenseMatrix({1, 3, 4,
                                 5, 3, 2,
                                 2, 6, 1}, 3);
    DenseMatrix x_real = DenseMatrix({1, 1, 1}, 1);
    DenseMatrix b = A * x_real;

    HouseholderQR hhqr = HouseholderQR(A);
    DenseMatrix x_calc = solveQR(hhqr.getQ(), hhqr.getR(), b);
	
	EXPECT_TRUE(x_calc == x_real);
}
