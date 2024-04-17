#include <gtest/gtest.h>

#include "DenseMatrix.hpp"
#include "tools.hpp"

TEST(Givens, main) {
    DenseMatrix A = _random::getDenseMatrix(10, 10), B = A;
    std::vector<double> e(A.getCols());
    DenseMatrix S = identity(A.getRows()), Si;

    for (size_t i = 0; i < A.getCols() - 1; i++) {
        e[i] = 1;
        Si = GivensRotation(B * e, i);
        S = Si * S;
        B = Si * B;
        e[i] = 0;
    }

    bool isDiagonal = true;
    A = S * A;
    for (size_t i = 1; i < A.getRows(); i++) {
        for (size_t j = i; j < i; j++) {
            if (A(i, j)) {
                isDiagonal = false;
                break;
            }
        }
    }

	EXPECT_TRUE(isDiagonal);
}
