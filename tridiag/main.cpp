#include <gtest/gtest.h>
#include "tridiag_solver.hpp"

TEST(tridiag, main) {
	std::vector<double> calc_solution = solve_tridiag("example01");
	std::vector<double> true_solution = get_solution("example01");
	
	ASSERT_NEAR(calc_solution[0], true_solution[0], 1e-16);
}