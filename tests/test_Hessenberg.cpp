#include <gtest/gtest.h>

#include "HessenbergMatrix.hpp"
#include "DenseMatrix.hpp"
#include "tools.hpp"

TEST(Hessenberg, main) {
    unsigned n = 5;
    // CSRMatrix A({1, 1}, {0, 1}, {0, 0, 0, 1, 2, 2}); 
    CSRMatrix A = _random::getCSRMatrix(n, n, 0.1);
    std::vector<double> x(n, 1);
    std::vector<double> b = A * x;
    std::vector<double> x0(n);

    HessenbergMatrix hm(A, b, x0);

    unsigned m = 4;
    for (unsigned i = 0; i < m; i++) {
        hm.iterate();
        std::cout << hm.getResidual() << std::endl;
    }

    /*
    std::cout << A << std::endl;
    std::cout << hm.getR() << std::endl;
    std::cout << hm.getQ() << std::endl;
    std::cout << hm.getQ() * hm.getR() << std::endl;
    */
	EXPECT_TRUE(hm.getQ() * hm.getQ().transpose() == identity(m + 1));
}
