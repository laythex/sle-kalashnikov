#include "CSR.hpp"

CSR::CSR() {}

void CSR::initialize(const std::vector<double>& vals, const std::vector<unsigned int>& cols, const std::vector<unsigned int>& rows) {
    this->vals = vals;
    this->cols = cols;
    this->rows = rows;
}

void CSR::initialize(const std::vector<double>& matrix, unsigned int number_of_columns) {
    for (unsigned int i = 0; i < matrix.size(); i++) {
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

double CSR::operator()(unsigned int i, unsigned int j) const {
    for (unsigned int k = rows[i]; k < rows[i + 1]; k++) {
        if (cols[k] == j) {
            return vals[k];
        }
    }
    return 0;
}

std::vector<double> CSR::operator*(const std::vector<double>& v) const {
    std::vector<double> res(v.size());
    for (unsigned int i = 0; i < v.size(); i++) {
        for (unsigned int k = rows[i]; k < rows[i + 1]; k++) {
            res[i] += vals[k] * v[cols[k]];
        }
    }
    return res;
}

double CSR::operator*(const CSR& other) const {
    double res = 0;
    for (unsigned int i = 1; i < rows.size(); i++) {
        if ((rows[i] - rows[i - 1]) * (other.rows[i] - other.rows[i - 1])) {
            res += vals[rows[i - 1]] * other.vals[other.rows[i - 1]];
        }
    }
    return res;
}

CSR CSR::operator+(const CSR& other) const {
    CSR res;
    res.rows.push_back(0);

    for (unsigned int i = 1; i < rows.size(); i++) {
        unsigned int left = rows[i] - rows[i - 1];
        unsigned int right = other.rows[i] - other.rows[i - 1];

        if (left || right) {
            res.vals.push_back(left * vals[rows[i - 1]] + right * other.vals[other.rows[i - 1]]);
            res.cols.push_back(0);
        }

        res.rows.push_back(res.rows.back() + (left || right));
    }

    return res;
}

std::ostream& operator<<(std::ostream& os, const CSR& csr) {
    os << "Vals:\t{ "; 
    for (unsigned int i = 0; i < csr.vals.size(); i++) {
        os << csr.vals[i] << ' ';
    }
    os << "}" << std::endl;

    os << "Cols:\t{ "; 
    for (unsigned int i = 0; i < csr.cols.size(); i++) {
        os << csr.cols[i] << ' ';
    }
    os << "}" << std::endl;

    os << "Rows:\t{ "; 
    for (unsigned int i = 0; i < csr.rows.size(); i++) {
        os << csr.rows[i] << ' ';
    }
    os << "}" << std::endl;

    return os;
}