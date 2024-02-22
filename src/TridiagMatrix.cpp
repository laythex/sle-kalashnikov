#include "TridiagMatrix.hpp"

void TridiagMatrix::initialize(const std::vector<double>& a_short_data, const std::vector<double>& b_data, const std::vector<double>& c_short_data) {
	a = a_short_data;
	b = b_data;
	c = c_short_data;

	a.insert(a.begin(), 0);
	c.push_back(0);

	n = b.size();
}

// double TridiagMatrix::operator()(int i, int j) const {}
// std::vector<double> TridiagMatrix::voperator*(const std::vector<double>& v) const {}
// TridiagMatrix TridiagMatrix::operator+(const TridiagMatrix& other) const {}
// double TridiagMatrix::operator*(const TridiagMatrix& other) const {}
// std::ostream& operator<<(std::ostream& os, const TridiagMatrix& tm) {}

std::vector<double> TridiagMatrix::solve(const std::vector<double>& d) const {
	std::vector<double> p(n), q(n), x(n);

	for (int i = 0; i < n - 1; i++) {
		p[i + 1] = -c[i] / (a[i] * p[i] + b[i]);
		q[i + 1] = (d[i] - a[i] * q[i]) / (a[i] * p[i] + b[i]);
	}

	x[n - 1] = (d[n - 1] - a[n - 1] * q[n - 1]) / (a[n - 1] * p[n - 1] + b[n - 1]);

	for (int i = n - 2; i >= 0; i--) {
		x[i] = p[i + 1] * x[i + 1] + q[i + 1];
	}

	return x;
}
