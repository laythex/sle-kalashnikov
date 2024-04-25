#include <gtest/gtest.h>

#include "HessenbergMatrix.hpp"
#include "DenseMatrix.hpp"
#include "tools.hpp"

TEST(Hessenberg, main) {
    unsigned n = 10;
    CSRMatrix A = _random::getCSRMatrix(n, n, 0.1);
    std::vector<double> x = _random::getVector(n);
    std::vector<double> b = A * x;
    std::vector<double> x0(n), r0 = A * x0 - b;

    HessenbergMatrix hm(A, r0);
    hm.iterate();
    hm.iterate();
    hm.iterate();
    hm.iterate();

    bool isDiagonal = true;
    for (size_t i = 1; i < hm.H.getRows(); i++) {
        for (size_t j = i; j < i; j++) {
            if (hm.H(i, j)) {
                isDiagonal = false;
                break;
            }
        }
    }
    
	EXPECT_TRUE(isDiagonal);
}
