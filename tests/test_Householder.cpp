#include <gtest/gtest.h>

#include "DenseMatrix.hpp"
#include "tools.hpp"
#include "HouseholderQR.hpp"

TEST(Householder, main) {
	DenseMatrix A = _random::getDenseMatrix(100, 100, 1, -1, 1);
    HouseholderQR hhqr = HouseholderQR(A);
    DenseMatrix Q = hhqr.getQ();
    DenseMatrix R = hhqr.getR();
	
	EXPECT_TRUE(Q * R == A);
}
