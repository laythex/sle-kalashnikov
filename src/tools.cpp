#include "tools.hpp"

double random_double(double min, double max) {
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    static std::default_random_engine e(seed);
    std::uniform_real_distribution<double> d(min, max);
    return d(e);
}

DenseMatrix randomDenseMatrix(unsigned rows, unsigned cols, double density, double min_el, double max_el, bool integer) {
    std::vector<double> data(rows * cols);
    double el;
    for (unsigned i = 0; i < rows * cols; i++) {
        if (random_double(0, 1) < density) {
            el = random_double(min_el, max_el);
            data[i] = integer ? trunc(el) : el;
        }
    }
    return DenseMatrix(data, cols);
}

CSRMatrix randomCSRMatrix(unsigned rows, unsigned cols, double density, double min_el, double max_el) {
    std::vector<double> data(rows * cols);
    for (unsigned i = 0; i < rows * cols; i++) {
        if (random_double(0, 1) < density) {
            data[i] = random_double(min_el, max_el);
        }
    }

    CSRMatrix csr;
    csr.initialize(data, cols);
    return csr;
}

std::vector<double> randomVector(unsigned rows, double min_el, double max_el) {
    std::vector<double> res(rows);
    for (unsigned i = 0; i < rows; i++) {
        res[i] = random_double(min_el, max_el);
    }
    return res;
}
