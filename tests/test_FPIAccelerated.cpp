#include <gtest/gtest.h>

#include "CSRMatrix.hpp"
#include "tools.hpp"
#include "solvers.hpp"

TEST(FPIAccelerated, main) {
    CSRMatrix A = CSRMatrix({2, 1, 5, 4, 3}, {0, 1, 1, 0, 2}, {0, 2, 3, 5});
	std::vector<double> x_real = {1, 1, 1}, x_calc;
    std::vector<double> b = A * x_real;

    unsigned m = 128;
    std::vector<double> roots = FPIAcceleratedTools::calcChebyshevRoots(m);
    std::vector<size_t> permutations = FPIAcceleratedTools::calcPermutations(m); 
    double lambda = FPIAcceleratedTools::calcMaxEigenvalue(A, 1e-3);

    x_calc = solvers::FPIAccelerated(A, b, lambda, roots, permutations, {0, 0, 0}, 1e-12);

    EXPECT_TRUE(x_calc == x_real);
}
