#include "HouseholderQR.hpp"

HouseholderQR::HouseholderQR(DenseMatrix dm) : A(dm) {
    factorize();
} 

DenseMatrix& HouseholderQR::getQ() {
    return Q;
}

DenseMatrix& HouseholderQR::getR() {
    return R;
}

std::vector<double> HouseholderQR::theta(const std::vector<double>& x, const std::vector<double>& v) const {
    return x - v * static_cast<double>(v * x) / (v * v) * 2.0;  
}

void HouseholderQR::factorize() {
    unsigned n = A.getRows();
    
    std::vector<double> rawR = A.transpose().getData(),
                        rawQ = identity(n).getData();

    for (unsigned i = 0; i < n - 1; i++) {
        std::vector<double> x = std::vector<double>(rawR.begin() + n * i + i, rawR.begin() + n * (i + 1));
        std::vector<double> e(n - i);
        e[0] = 1;
        double sign = x[i] >= 0 ? 1 : -1;
        std::vector<double> v = x + e * norm2(x) * sign;

        std::vector<double> tmp(n - i);
        for (unsigned j = i; j < n; j++) {
            tmp = std::vector<double>(rawR.begin() + n * j + i, rawR.begin() + n * (j + 1));
            tmp = theta(tmp, v);
            for (unsigned k = n * j + i; k < n * (j + 1); k++) {
                rawR[k] = tmp[k - (n * j + i)];
            }
        }

        for (unsigned j = 0; j < n; j++) { 
            tmp = std::vector<double>(rawQ.begin() + n * j + i, rawQ.begin() + n * (j + 1));
            tmp = theta(tmp, v);
            for (unsigned k = n * j + i; k < n * (j + 1); k++) {
                rawQ[k] = tmp[k - (n * j + i)];
            }
        } 
    }

    R = DenseMatrix(rawR, n).transpose();
    Q = DenseMatrix(rawQ, n);
}
