#include "tools.hpp"

std::vector<double> operator*(const std::vector<double>& left, double right) {
    unsigned n = left.size();
    std::vector<double> res(n);

    for (unsigned i = 0; i < left.size(); i++) {
        res[i] = left[i] * right;
    }

    return res;
}

std::vector<double> operator*(double left, const std::vector<double>&  right) {
    return right * left;
}

std::vector<double> operator+(const std::vector<double>& left, const std::vector<double>& right) {
    unsigned n = left.size();
    std::vector<double> res(n);

    for (unsigned i = 0; i < left.size(); i++) {
        res[i] = left[i] + right[i];
    }

    return res;
}

std::vector<double> operator-(const std::vector<double>& left, const std::vector<double>& right) {
    return left + (-1) * right;
}

double operator*(const std::vector<double>& left, const std::vector<double>& right) {
    unsigned n = left.size();
    double res = 0;

    for (unsigned i = 0; i < left.size(); i++) {
        res += left[i] * right[i];
    }

    return res;
}

bool operator==(const std::vector<double>& left, const std::vector<double>& right) {
    double eps = 1e-13;

    for (unsigned i = 0; i < left.size(); i++) {
        if (fabs(left[i] - right[i]) > eps) {
            return false;
        } 
    }

    return true;
}

double norm2(const std::vector<double>& vec) {
    return sqrt(vec * vec);
}

std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec) {
    for (unsigned i = 0; i < vec.size(); i++) {
        os << vec[i] << ' ';
    }

    return os;
}

double _random::getDouble(double min, double max) {
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    static std::default_random_engine e(seed);
    std::uniform_real_distribution<double> d(min, max);
    return d(e);
}

DenseMatrix _random::getDenseMatrix(unsigned rows, unsigned cols, double density, double min_el, double max_el, bool integer) {
    std::vector<double> data(rows * cols);
    double el;

    for (unsigned i = 0; i < rows * cols; i++) {
        if (_random::getDouble(0, 1) < density) {
            el = _random::getDouble(min_el, max_el);
            data[i] = integer ? trunc(el) : el;
        }
    }

    return DenseMatrix(data, cols);
}

CSRMatrix _random::getCSRMatrix(unsigned rows, unsigned cols, double density, double min_el, double max_el) {
    std::vector<double> data(rows * cols);

    for (unsigned i = 0; i < rows * cols; i++) {
        if (_random::getDouble(0, 1) < density) {
            data[i] = _random::getDouble(min_el, max_el);
        }
    }

    return CSRMatrix(data, cols);
}

std::vector<double> _random::getVector(unsigned rows, double min_el, double max_el) {
    std::vector<double> res(rows);

    for (unsigned i = 0; i < rows; i++) {
        res[i] = _random::getDouble(min_el, max_el);
    }

    return res;
}

std::vector<double> Jacobi::multiply(const CSRMatrix& csr, const std::vector<double>& v) {
    std::vector<double> res(v.size());

    for (unsigned i = 0; i < v.size(); i++) {
        for (unsigned j = csr.rowsAt(i); j < csr.rowsAt(i + 1); j++) {
            if (i != csr.colsAt(j)) {
                res[i] += csr.valsAt(j) * v[csr.colsAt(j)];
            }
        }
    }

    return res;
}

CSRMatrix Jacobi::inverseDiagonal(const CSRMatrix& csr) {
    unsigned n = csr.getRowsSize() - 1;
    std::vector<double> vals(n);
    std::vector<unsigned> cols(n), rows(n + 1);

    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = csr.rowsAt(i); j < csr.rowsAt(i + 1); j++) {
            if (i == csr.colsAt(j)) {
                vals[i] = 1 / csr.valsAt(j);
                cols[i] = i;
                rows[i + 1] = i + 1;
            }
        }
    }

    return CSRMatrix(vals, cols, rows);
}

std::vector<double> GaussSeidel::inverseDiagonal(const CSRMatrix& csr) {
    unsigned n = csr.getRowsSize() - 1;
    std::vector<double> res(n);

    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = csr.rowsAt(i); j < csr.rowsAt(i + 1); j++) {
            if (i == csr.colsAt(j)) {
                res[i] = 1 / csr.valsAt(j);
            }
        }
    }

    return res;
}
