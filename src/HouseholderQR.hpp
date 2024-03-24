#pragma once

#include "DenseMatrix.hpp"
#include "tools.hpp"

class HouseholderQR {
    public:
        HouseholderQR(DenseMatrix dm); 
        
        DenseMatrix& getQ();
        DenseMatrix& getR();

    protected:
        void factorize();
        std::vector<double> theta(const std::vector<double>& x, const std::vector<double>& v) const;

        DenseMatrix A, Q, R;
};