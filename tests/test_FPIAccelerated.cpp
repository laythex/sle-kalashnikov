#include <gtest/gtest.h>

#include "CSRMatrix.hpp"
#include "tools.hpp"
#include "solvers.hpp"

TEST(FPIAccelerated, main) {
    CSRMatrix A = CSRMatrix({2, 1, 5, 4, 3}, {0, 1, 1, 0, 2}, {0, 2, 3, 5});
	std::vector<double> x_real = {1, 1, 1}, x_calc;
    std::vector<double> b = A * x_real;

    x_calc = solvers::FPIAccelerated(A, b, 128, {0, 0, 0}, 1e-12);

    EXPECT_TRUE(x_calc == x_real);
}
