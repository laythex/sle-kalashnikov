#include "DenseMatrix.hpp"

DenseMatrix::DenseMatrix() {}

void DenseMatrix::initialize(const std::vector<double>& matrix_data, unsigned number_of_columns) {
    data = matrix_data;
    cols = number_of_columns;
    rows = data.size() / cols;
}

double DenseMatrix::operator()(int i, int j) const {
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
    std::vector<double> res_matrix(data.size());
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            res_matrix[i * cols + j] += (*this)(i, j) + other(i, j);
        }
    }
    DenseMatrix res;
    res.initialize(res_matrix, cols);
    return res;
}

double DenseMatrix::operator*(const DenseMatrix& other) const {
    double res = 0;
    for (unsigned i = 0; i < rows; i++) {
        res += (*this)(i, 0) * other(i, 0);
    }
    return res;
}

std::ostream& operator<<(std::ostream& os, const DenseMatrix& dm) {
    for (unsigned i = 0; i < dm.rows; i++) {
        for (unsigned j = 0; j < dm.cols; j++) {
            os << dm(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}
