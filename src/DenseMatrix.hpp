#pragma once

#include <cmath>
#include <vector>
#include <iostream>

class DenseMatrix {
    public:
        DenseMatrix();
        DenseMatrix(const std::vector<double>& data, size_t cols);
        
        double operator()(size_t i, size_t j) const;
        double& at(size_t i, size_t j);
        std::vector<double> operator*(const std::vector<double>& v) const;
        DenseMatrix operator+(const DenseMatrix& other) const;
        DenseMatrix operator-(const DenseMatrix& other) const;
        DenseMatrix operator*(const DenseMatrix& other) const;
        DenseMatrix operator*(double x) const;
        DenseMatrix operator/(double x) const;
        bool operator==(const DenseMatrix& other) const;
        DenseMatrix transpose() const;
        DenseMatrix resize(size_t nrows, size_t ncols) const;

        size_t getRows() const;
        size_t getCols() const;
        const std::vector<double>& getData() const;
        const std::vector<double> getRow(size_t k) const;
        const std::vector<double> getCol(size_t k) const;
    
    private:
        size_t rows, cols;
        std::vector<double> data;
};

std::ostream& operator<<(std::ostream& os, const DenseMatrix& dm);

DenseMatrix identity(size_t n, size_t m);
DenseMatrix identity(size_t n);

DenseMatrix zeros(size_t n, size_t m);
DenseMatrix zeros(size_t n);

DenseMatrix ones(size_t n, size_t m);
DenseMatrix ones(size_t n);