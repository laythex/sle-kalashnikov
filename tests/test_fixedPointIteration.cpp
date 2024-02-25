#include <gtest/gtest.h>

#include "DenseMatrix.hpp"
#include "tools.hpp"
#include "solvers.hpp"

TEST(fixedPointIteration, main) {
    CSRMatrix A = CSRMatrix({2, 1, 5, 4, 3}, {0, 1, 1, 0, 2}, {0, 2, 3, 5});
	std::vector<double> x_real = {1, 1, 1}, x_calc;
    std::vector<double> b = A * x_real;

    x_calc = solvers::fixedPointIteration(A, b, 0.1, {-15, 1e3, 1e-5}, 1e-13);

    std::cout << x_calc;
    EXPECT_TRUE(x_calc == x_real);
}
