#include "CSRMatrix.hpp"

CSRMatrix::CSRMatrix() : vals({}), cols({}), rows({}) {}

CSRMatrix::CSRMatrix(const std::vector<double>& vals, const std::vector<unsigned>& cols, const std::vector<unsigned>& rows) : vals(vals), cols(cols), rows(rows) {}

CSRMatrix::CSRMatrix(const std::vector<double>& matrix, unsigned number_of_columns) {
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

double CSRMatrix::valsAt(unsigned i) const {
    return vals[i];
}

unsigned CSRMatrix::colsAt(unsigned i) const {
    return cols[i];
}

unsigned CSRMatrix::rowsAt(unsigned i) const {
    return rows[i];
}

unsigned CSRMatrix::getValsSize() const {
    return vals.size();
}

unsigned CSRMatrix::getColsSize() const {
    return cols.size();
}

unsigned CSRMatrix::getRowsSize() const {
    return rows.size();
}

std::ostream& operator<<(std::ostream& os, const CSRMatrix& csr) {
    os << "Vals = { "; 
    for (unsigned i = 0; i < csr.getValsSize(); i++) {
        os << csr.valsAt(i) << ' ';
    }
    os << "}" << '\n';

    os << "Cols = { "; 
    for (unsigned i = 0; i < csr.getColsSize(); i++) {
        os << csr.colsAt(i) << ' ';
    }
    os << "}" << '\n';

    os << "Rows = { "; 
    for (unsigned i = 0; i < csr.getRowsSize(); i++) {
        os << csr.rowsAt(i) << ' ';
    }
    os << "}" << '\n';

    return os;
}