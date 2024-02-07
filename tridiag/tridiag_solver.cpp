#include <fstream>
#include <vector>
#include <gtest/gtest.h>

std::vector<double> solve_tridiag(std::string filename) {
	std::ifstream file;
	file.open("../examples/tridiag/" + filename + "/data.txt");

	int n;
	file >> n;
	
	std::vector<double> matrix(n * 7, 0);
	for (int i = 0; i < n * 4; i++) {
		if (i != 0 && i != n * 3 - 1) { // a_0 and c_n
			file >> matrix[i];
		}
	}

	for (int i = 0; i < n - 1; i++) {
		matrix[n * 4 + i + 1] = -matrix[n * 2 + i] / (matrix[i] * matrix[n * 4 + i] + matrix[n + i]);
		matrix[n * 5 + i + 1] = (matrix[n * 3 + i] - matrix[i] * matrix[n * 5 + 1]) / (matrix[i] * matrix[n * 4 + i] + matrix[n + i]);
	}

	matrix[n * 7 - 1] = (matrix[n * 4 - 1] - matrix[n - 1] * matrix[n * 6 - 1]) / (matrix[n - 1] * matrix[n * 5 - 1] + matrix[n * 2 - 1]);

	for (int i = n - 2; i >= 0; i--) {
		matrix[n * 6 + i] = matrix[n * 4 + i + 1] * matrix[n * 6 + i + 1] + matrix[n * 5 + i + 1];
	}

	return std::vector<double>(matrix.begin() + n * 6, matrix.end());
}

std::vector<double> get_solution(std::string filename) {
	std::ifstream file;
	file.open("../examples/tridiag/" + filename + "/solution.txt");

	std::vector<double> solution;
	double element;

	while (file >> element) {
		solution.push_back(element);
	}

	for (int i = 0; i < 3; i++) {
		std::cout << solution[i] << std::endl;
	}

	return solution;
}

TEST(tridiag, main) {
	std::vector<double> calc_solution = solve_tridiag("example01");
	std::vector<double> true_solution = get_solution("example01");
	
	ASSERT_NEAR(calc_solution[0], true_solution[0], 1e-16);
}