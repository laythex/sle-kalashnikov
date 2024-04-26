#include "HessenbergMatrix.hpp"

HessenbergMatrix::HessenbergMatrix(const CSRMatrix& A, const std::vector<double>& b, const std::vector<double>& x0) : n(A.getRowsSize()), j(0), A(A), b(b), x0(x0) {
    std::vector<double> r0 = A * x0 - b;
    r0_norm = norm2(r0);
    V.push_back(r0 / r0_norm); // Находим первый базисный вектор
    q.push_back(1); // Начинаем строить первый столбец Q^T
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

    // Применяем все предыдущие повороты к новому столбцу матрицы H
    double h1, h2;
    for (size_t i = 0; i < j; i++) {
        h1 = rots[i].c * H(i, j) - rots[i].s * H(i + 1, j);
        h2 = rots[i].s * H(i, j) + rots[i].c * H(i + 1, j);
        H.at(i, j) = h1;
        H.at(i + 1, j) = h2;
    }

    // Вычисляем новый поворот
    double tmp = sqrt(H(j, j) * H(j, j) + h * h);
    double c = H(j, j) / tmp;
    double s = -h / tmp;
    rots.push_back({c, s}); // Запоминаем поворот

    // Превращаем матрицу H в матрицу R
    H.at(j, j) = tmp;

    // Вычисляем первый столбец Q^T
    q.push_back(0);
    h1 = c * q[j] - s * q[j + 1];
    h2 = s * q[j] + c * q[j + 1];
    q[j] = h1;
    q[j + 1] = h2;

    j++;
}

double HessenbergMatrix::getResidual() const {
    return fabs(q[j]) * r0_norm;
}

std::vector<double> HessenbergMatrix::solve() const {
    std::vector<double> dx(x0.size());
    std::vector<double> y(j);
    std::vector<double> z = q * r0_norm;

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
