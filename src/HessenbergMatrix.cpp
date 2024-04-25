#include "HessenbergMatrix.hpp"

HessenbergMatrix::HessenbergMatrix(const CSRMatrix& A, const std::vector<double>& r0) : n(A.getRowsSize()), j(0), A(A), r0(r0) { 
    V.push_back(normalize(r0));
}

void HessenbergMatrix::iterate() {
    std::vector<double> t = A * V[j];
    H = H.resize(j + 2, j + 1);

    double h;
    for (size_t k = 0; k < j + 1; k++) {
        h = V[k] * t;
        t = t - h * V[k];
        H.at(k, j) = h; 
    }

    h = norm2(t);
    V.push_back(t / h);

    rots.push_back(GivensRotation(H.getCol(j), j));
    H.at(j, j) = sqrt(H(j, j) * H(j, j) + h * h);

    j++;
}