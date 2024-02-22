#pragma once

#include "DenseMatrix.hpp"

class HouseholderQR {
    public:
        HouseholderQR(DenseMatrix dm); 
        
        DenseMatrix getQ() const;
        DenseMatrix getR() const;

    protected:
        void factorize();
        DenseMatrix theta(const DenseMatrix& x, const DenseMatrix& v) const;

        DenseMatrix A, Q, R;
};