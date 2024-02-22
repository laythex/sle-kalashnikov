#include "CSRMatrix.hpp"

void CSRMatrix::initialize(const std::vector<double>& vals, const std::vector<unsigned>& cols, const std::vector<unsigned>& rows) {
    this->vals = vals;
    this->cols = cols;
    this->rows = rows;
}

void CSRMatrix::initialize(const std::vector<double>& matrix, unsigned number_of_columns) {
    for (unsigned i = 0; i < matrix.size(); i++) {
        if (i % number_of_columns == 0) {
            rows.push_back(vals.size());
        }

        if (matrix[i]) {
            vals.push_back(matrix[i]);
            cols.push_back(i % number_of_columns);
        }
    }
    rows.push_back(vals.size());
}

double CSRMatrix::operator()(unsigned i, unsigned j) const {
    for (unsigned k = rows[i]; k < rows[i + 1]; k++) {
        if (cols[k] == j) {
            return vals[k];
        }
    }
    return 0;
}

std::vector<double> CSRMatrix::operator*(const std::vector<double>& v) const {
    std::vector<double> res(v.size());
    for (unsigned i = 0; i < v.size(); i++) {
        for (unsigned k = rows[i]; k < rows[i + 1]; k++) {
            res[i] += vals[k] * v[cols[k]];
        }
    }
    return res;
}

double CSRMatrix::operator*(const CSRMatrix& other) const {
    double res = 0;
    for (unsigned i = 1; i < rows.size(); i++) {
        if ((rows[i] - rows[i - 1]) * (other.rows[i] - other.rows[i - 1])) {
            res += vals[rows[i - 1]] * other.vals[other.rows[i - 1]];
        }
    }
    return res;
}

CSRMatrix CSRMatrix::operator+(const CSRMatrix& other) const {
    CSRMatrix res;
    res.rows.push_back(0);

    for (unsigned i = 1; i < rows.size(); i++) {
        unsigned left = rows[i] - rows[i - 1];
        unsigned right = other.rows[i] - other.rows[i - 1];

        if (left || right) {
            res.vals.push_back(left * vals[rows[i - 1]] + right * other.vals[other.rows[i - 1]]);
            res.cols.push_back(0);
        }

        res.rows.push_back(res.rows.back() + (left || right));
    }

    return res;
}

std::ostream& operator<<(std::ostream& os, const CSRMatrix& csr) {
    os << "Vals = { "; 
    for (unsigned i = 0; i < csr.vals.size(); i++) {
        os << csr.vals[i] << ' ';
    }
    os << "}" << std::endl;

    os << "Cols = { "; 
    for (unsigned i = 0; i < csr.cols.size(); i++) {
        os << csr.cols[i] << ' ';
    }
    os << "}" << std::endl;

    os << "Rows = { "; 
    for (unsigned i = 0; i < csr.rows.size(); i++) {
        os << csr.rows[i] << ' ';
    }
    os << "}" << std::endl;

    return os;
}