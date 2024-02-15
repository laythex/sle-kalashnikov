#pragma once

#include <vector>
#include <iostream>

class CSR {
    public:
        CSR();

        void initialize(const std::vector<double>& vals, const std::vector<unsigned int>& cols, const std::vector<unsigned int>& rows);
        void initialize(const std::vector<double>& matrix, unsigned int number_of_columns);

        double operator()(unsigned int i, unsigned int j) const;
        std::vector<double> operator*(const std::vector<double>& v) const;
        CSR operator+(const CSR& other) const;
        double operator*(const CSR& other) const;
        friend std::ostream& operator<<(std::ostream& os, const CSR& csr);
        
    protected:
        std::vector<double> vals;
        std::vector<unsigned int> cols;
        std::vector<unsigned int> rows;
};


