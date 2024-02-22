#include "HouseholderQR.hpp"

HouseholderQR::HouseholderQR(DenseMatrix dm) : A(dm) {
    factorize();
} 

DenseMatrix HouseholderQR::getQ() const {
    return Q;
}

DenseMatrix HouseholderQR::getR() const {
    return R;
}

DenseMatrix HouseholderQR::theta(const DenseMatrix& x, const DenseMatrix& v) const {
    return x - v * static_cast<double>(v.transpose() * x) / (v.transpose() * v) * 2.0;  
}

void HouseholderQR::factorize() {
    unsigned n = A.getRows();
    
    std::vector<double> rawR = A.transpose().getData(),
                        rawQ = identity(n).getData(), 
                        tmp;

    for (unsigned i = 0; i < n - 1; i++) {
        tmp = std::vector<double>(rawR.begin() + n * i + i, rawR.begin() + n * (i + 1));
        DenseMatrix x = DenseMatrix(tmp, 1);

        DenseMatrix e = identity(n - i, 1);
        double sign = x(i, 0) >= 0 ? 1 : -1;
        DenseMatrix v = x + e * x.norm() * sign;

        DenseMatrix res;
        for (unsigned j = 0; j < n; j++) {
            tmp = std::vector<double>(rawR.begin() + n * j + i, rawR.begin() + n * (j + 1));
            res = theta(DenseMatrix(tmp, 1), v);
            for (unsigned k = n * j + i; k < n * (j + 1); k++) {
                rawR[k] = res(k - (n * j + i), 0);
            }

            tmp = std::vector<double>(rawQ.begin() + n * j + i, rawQ.begin() + n * (j + 1));
            res = theta(DenseMatrix(tmp, 1), v);
            for (unsigned k = n * j + i; k < n * (j + 1); k++) {
                rawQ[k] = res(k - (n * j + i), 0);
            }
        }
    }

    R = DenseMatrix(rawR, n).transpose();
    Q = DenseMatrix(rawQ, n);
}
