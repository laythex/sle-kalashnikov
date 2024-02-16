#include <vector>
#include <gtest/gtest.h>

#include "TridiagMatrix.hpp"

TEST(tridiag, main) {
	std::vector <double> a = {1, 3};
	std::vector <double> b = {4, 6, 7};
	std::vector <double> c = {1, 3};
	std::vector <double> d = {1, 1, 1};
	std::vector <double> x_real = {0.232, 0.072, 0.112};
	std::vector <double> x_calc;

	TridiagMatrix tm;
	tm.initialize(a, b, c);
	x_calc = tm.solve(d);

	ASSERT_NEAR(x_calc[0], x_real[0], 1e-16);
}
