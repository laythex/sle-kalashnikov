#pragma once

#include <vector>
#include <iostream>

class CSRMatrix {
    public:
        CSRMatrix();

        void initialize(const std::vector<double>& vals, const std::vector<unsigned>& cols, const std::vector<unsigned>& rows);
        void initialize(const std::vector<double>& matrix, unsigned number_of_columns);

        double operator()(unsigned i, unsigned j) const;
        std::vector<double> operator*(const std::vector<double>& v) const;
        CSRMatrix operator+(const CSRMatrix& other) const;
        double operator*(const CSRMatrix& other) const;
        friend std::ostream& operator<<(std::ostream& os, const CSRMatrix& csr);
        
    protected:
        std::vector<double> vals;
        std::vector<unsigned> cols;
        std::vector<unsigned> rows;
};


