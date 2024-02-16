#pragma once

#include <vector>
#include <iostream>

class DenseMatrix {
    public:
        DenseMatrix();

        void initialize(const std::vector<double>& matrix_data, unsigned number_of_columns);
        
        double operator()(int i, int j) const;
        std::vector<double> operator*(const std::vector<double>& v) const;
        DenseMatrix operator+(const DenseMatrix& other) const;
        double operator*(const DenseMatrix& other) const;
        friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& dm);
    
    protected:
        unsigned rows, cols;
        std::vector<double> data;
};