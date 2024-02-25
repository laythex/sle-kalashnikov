#include "DenseMatrix.hpp"

DenseMatrix::DenseMatrix() : rows(0),  cols(0), data({}) {}
DenseMatrix::DenseMatrix(const std::vector<double>& data, unsigned cols) : rows(data.size() / cols), cols(cols), data(data) {}

double DenseMatrix::operator()(unsigned i, unsigned j) const {
    return data[i * cols + j];
}

std::vector<double> DenseMatrix::operator*(const std::vector<double>& v) const {
    std::vector<double> res(v.size());
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            res[i] += (*this)(i, j) * v[j];
        }
    }
    return res;
}

DenseMatrix DenseMatrix::operator+(const DenseMatrix& other) const {
    std::vector<double> v(data.size());
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            v[i * cols + j] += (*this)(i, j) + other(i, j);
        }
    }
    return DenseMatrix(v, cols);
}

DenseMatrix DenseMatrix::operator-(const DenseMatrix& other) const {
    return (*this) + other * (-1.0);
}

DenseMatrix DenseMatrix::operator*(const DenseMatrix& other) const {
    unsigned n = rows, m = other.cols;
    std::vector<double> v(n * m);
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < m; j++) {
            for (unsigned k = 0; k < cols; k++) {
                v[m * i + j] += (*this)(i, k) * other(k, j);
            }
        }
    }
    return DenseMatrix(v, m);
}

DenseMatrix DenseMatrix::operator*(double x) const {
    std::vector<double> v(rows * cols);
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            v[i * cols + j] = (*this)(i, j) * x;
        }
    }
    return DenseMatrix(v, cols);
}

DenseMatrix DenseMatrix::operator/(double x) const {
    return (*this) * (1 / x);
}

bool DenseMatrix::operator==(const DenseMatrix& other) const {
    double eps = 1e-13;
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            if (fabs((*this)(i, j) - other(i, j)) > eps) {
                return false;
            } 
        }
    }
    return true;
}

DenseMatrix::operator double() const {
    return (*this)(0, 0);
}

DenseMatrix DenseMatrix::transpose() const {
    std::vector<double> v(rows * cols);
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            v[rows * j + i] = (*this)(i, j);
        }
    }
    return DenseMatrix(v, rows);
}

double DenseMatrix::norm() const {
    double norm = 0;
    for (unsigned i = 0; i < rows; i++) {
        norm += (*this)(i, 0) * (*this)(i, 0);
    }
    return std::sqrt(norm);
    
}

unsigned DenseMatrix::getRows() const {
    return rows;
}

unsigned DenseMatrix::getCols() const {
    return cols;
}

const std::vector<double>& DenseMatrix::getData() const {
    return data;
}

std::ostream& operator<<(std::ostream& os, const DenseMatrix& dm) {
    double eps = 1e-13, el;
    for (unsigned i = 0; i < dm.getRows(); i++) {
        for (unsigned j = 0; j < dm.getCols(); j++) {
            el = fabs(dm(i, j)) > eps ? dm(i, j) : 0;
            os << el << " ";
        }
        os << std::endl;
    }
    return os;
}    

DenseMatrix identity(unsigned n, unsigned m) {
    std::vector<double> v(n * m);
    for (unsigned i = 0; i < std::min(n, m); i++) {
        v[m * i + i] = 1;
    }
    return DenseMatrix(v, m);
}

DenseMatrix identity(unsigned n) {
    return identity(n, n);
}

DenseMatrix zeros(unsigned n, unsigned m) {
    std::vector<double> v(n * m);
    return DenseMatrix(v, m);
}

DenseMatrix zeros(unsigned n) {
    return zeros(n, n);
}

DenseMatrix ones(unsigned n, unsigned m) {
    std::vector<double> v(n * m, 1.0);
    return DenseMatrix(v, m);
}

DenseMatrix ones(unsigned n) {
    return ones(n, n);
}