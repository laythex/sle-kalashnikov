#pragma once

#include <vector>
#include <iostream>

class CSRMatrix {
    public:
        CSRMatrix();
        CSRMatrix(const std::vector<double>& vals, const std::vector<unsigned>& cols, const std::vector<unsigned>& rows);
        CSRMatrix(const std::vector<double>& matrix, unsigned number_of_columns);

        double operator()(unsigned i, unsigned j) const;
        std::vector<double> operator*(const std::vector<double>& v) const;

        double valsAt(unsigned i) const;
        unsigned colsAt(unsigned i) const;
        unsigned rowsAt(unsigned i) const;
    
        unsigned getValsSize() const;
        unsigned getColsSize() const;
        unsigned getRowsSize() const;
        
    private:
        std::vector<double> vals;
        std::vector<unsigned> cols;
        std::vector<unsigned> rows;
};

std::ostream& operator<<(std::ostream& os, const CSRMatrix& csr);
