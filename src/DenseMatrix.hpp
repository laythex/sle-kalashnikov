#pragma once

#include <cmath>
#include <vector>
#include <iostream>

class DenseMatrix {
    public:
        DenseMatrix();
        DenseMatrix(const std::vector<double>& data, unsigned cols);
        
        double operator()(unsigned i, unsigned j) const;
        std::vector<double> operator*(const std::vector<double>& v) const;
        DenseMatrix operator+(const DenseMatrix& other) const;
        DenseMatrix operator-(const DenseMatrix& other) const;
        DenseMatrix operator*(const DenseMatrix& other) const;
        DenseMatrix operator*(double x) const;
        DenseMatrix operator/(double x) const;
        bool operator==(const DenseMatrix& other) const;
        operator double() const;
        DenseMatrix transpose() const;
        double norm() const;

        unsigned getRows() const;
        unsigned getCols() const;
        const std::vector<double>& getData() const;
    
    private:
        unsigned rows, cols;
        std::vector<double> data;
};

std::ostream& operator<<(std::ostream& os, const DenseMatrix& dm);

DenseMatrix identity(unsigned n, unsigned m);
DenseMatrix identity(unsigned n);

DenseMatrix zeros(unsigned n, unsigned m);
DenseMatrix zeros(unsigned n);

DenseMatrix ones(unsigned n, unsigned m);
DenseMatrix ones(unsigned n);