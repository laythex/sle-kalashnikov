#include "HessenbergMatrix.hpp"

HessenbergMatrix::HessenbergMatrix(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0) : n(A.getRowsSize()), j(0), A(A), b(b), x0(x0) {
    std::vector<double> r0 = A * x0 - b;

    r0_norm = norm2(r0);
    V.push_back(r0 / r0_norm);

    Q = identity(1);
}

void HessenbergMatrix::iterate() {
    std::vector<double> t = A * V[j];
    H = H.resize(j + 2, j + 1);

    // Дополняем матрицу Хайзенберга
    double h;
    for (size_t k = 0; k < j + 1; k++) {
        h = V[k] * t;
        t = t - h * V[k];
        H.at(k, j) = h; 
    }
    h = norm2(t);
    V.push_back(t / h);

    // Считаем поворот T (косинус и синус)
    double eps = 1e-13;
    double tmp = sqrt(H(j, j) * H(j, j) + h * h);
    double c = tmp > eps ? H(j, j) / tmp : 1;       // Eсли столбец получился нулевым,
    double s = tmp > eps ? -h / tmp : 0;            // то не вращаем?

    std::cout << H(j, j) << ' ' << h << ' ' << tmp << std::endl;
    H.at(j, j) = tmp; // Матрицу H сразу превращаем в R

    // Считаем Q_{j} = T * Q_{j - 1}
    Q = Q.resize(j + 2, j + 2);
    Q.at(j + 1, j + 1) = 1;

    std::vector<double> q1 = Q.getRow(j); // Постарался максимально упростить матричное умножение
    for (size_t i = 0; i < j + 2; i++) {
        Q.at(j, i) = c * q1[i] - s * (i == j + 1);
        Q.at(j + 1, i) = s * q1[i] + c * (i == j + 1);
    }
    std::cout << "Q: " << std::endl;
    std::cout << getQ() << std::endl;
    std::cout << "R: " << std::endl;
    std::cout << getR() << std::endl;
    std::cout << "QR: " << std::endl;
    std::cout << getQ() * getR() << std::endl;

    j++;
}

double HessenbergMatrix::getResidual() const {
    return fabs(Q(j, 0)) * r0_norm;
}

std::vector<double> HessenbergMatrix::solve() const {
    std::vector<double> dx(x0.size());
    std::vector<double> y(j);
    std::vector<double> z = Q.getCol(0) * r0_norm;

    double tmp;
    for (int i = j - 1; i >= 0; i--) {
        tmp = 0;
        for (size_t k = j - 1; k > i; k--) {
            tmp += H(i, k) * y[k];
        }

        y[i] = (z[i] - tmp) / H(i, i);
        dx = dx + V[i] * y[i];
    }

    return x0 - dx;
}

DenseMatrix HessenbergMatrix::getQ() const {
    return Q.transpose();
}

DenseMatrix HessenbergMatrix::getR() const {
    return H;
}
