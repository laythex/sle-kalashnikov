#include "DenseMatrix.hpp"

DenseMatrix::DenseMatrix() : rows(0),  cols(0), data({}) {}
DenseMatrix::DenseMatrix(const std::vector<double>& data, size_t cols) : rows(data.size() / cols), cols(cols), data(data) {}

double DenseMatrix::operator()(size_t i, size_t j) const {
    return data[i * cols + j];
}

double& DenseMatrix::at(size_t i, size_t j) {
    return data[i * cols + j];
}

std::vector<double> DenseMatrix::operator*(const std::vector<double>& v) const {
    std::vector<double> res(v.size());
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            res[i] += operator()(i, j) * v[j];
        }
    }
    return res;
}

DenseMatrix DenseMatrix::operator+(const DenseMatrix& other) const {
    std::vector<double> v(data.size());
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            v[i * cols + j] += operator()(i, j) + other(i, j);
        }
    }
    return DenseMatrix(v, cols);
}

DenseMatrix DenseMatrix::operator-(const DenseMatrix& other) const {
    return (*this) + other * (-1.0);
}

DenseMatrix DenseMatrix::operator*(const DenseMatrix& other) const {
    size_t n = rows, m = other.cols;
    std::vector<double> v(n * m);
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < cols; k++) {
            for (size_t j = 0; j < m; j++) {
                v[m * i + j] += operator()(i, k) * other(k, j);
            }
        }
    }
    return DenseMatrix(v, m);
}

DenseMatrix DenseMatrix::operator*(double x) const {
    std::vector<double> v(rows * cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            v[i * cols + j] = operator()(i, j) * x;
        }
    }
    return DenseMatrix(v, cols);
}

DenseMatrix DenseMatrix::operator/(double x) const {
    return (*this) * (1 / x);
}

bool DenseMatrix::operator==(const DenseMatrix& other) const {
    double eps = 1e-13;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            if (fabs(operator()(i, j) - other(i, j)) > eps) {
                return false;
            } 
        }
    }
    return true;
}

DenseMatrix DenseMatrix::transpose() const {
    std::vector<double> v(rows * cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            v[j * rows + i] = operator()(i, j);
        }
    }
    return DenseMatrix(v, rows);
}

DenseMatrix DenseMatrix::resize(size_t nrows, size_t ncols) const {
    std::vector<double> v(nrows * ncols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            v[i * ncols + j] = operator()(i, j);
        }
    }
    return DenseMatrix(v, ncols);
}

size_t DenseMatrix::getRows() const {
    return rows;
}

size_t DenseMatrix::getCols() const {
    return cols;
}

const std::vector<double>& DenseMatrix::getData() const {
    return data;
}

const std::vector<double> DenseMatrix::getRow(size_t k) const {
    std::vector<double> v(cols);

    for (size_t i = 0; i < cols; i++) {
        v[i] = operator()(k, i);
    }

    return v;
}

const std::vector<double> DenseMatrix::getCol(size_t k) const {
    std::vector<double> v(rows);
    
    for (size_t i = 0; i < rows; i++) {
        v[i] = operator()(i, k);
    }

    return v;
}

std::ostream& operator<<(std::ostream& os, const DenseMatrix& dm) {
    double eps = 1e-13, el;
    for (size_t i = 0; i < dm.getRows(); i++) {
        for (size_t j = 0; j < dm.getCols(); j++) {
            el = fabs(dm(i, j)) > eps ? dm(i, j) : 0;
            os << el << " ";
        }
        os << std::endl;
    }
    return os;
}    

DenseMatrix identity(size_t n, size_t m) {
    std::vector<double> v(n * m);
    for (size_t i = 0; i < std::min(n, m); i++) {
        v[m * i + i] = 1;
    }
    return DenseMatrix(v, m);
}

DenseMatrix identity(size_t n) {
    return identity(n, n);
}

DenseMatrix zeros(size_t n, size_t m) {
    std::vector<double> v(n * m);
    return DenseMatrix(v, m);
}

DenseMatrix zeros(size_t n) {
    return zeros(n, n);
}

DenseMatrix ones(size_t n, size_t m) {
    std::vector<double> v(n * m, 1.0);
    return DenseMatrix(v, m);
}

DenseMatrix ones(size_t n) {
    return ones(n, n);
}