#pragma once

#include <vector>
#include <array>
#include <chrono>
#include <random>

#include "DenseMatrix.hpp"
#include "CSRMatrix.hpp"

std::vector<double> operator*(const std::vector<double>& left, double right);
std::vector<double> operator*(double left, const std::vector<double>&  right);
std::vector<double> operator/(const std::vector<double>& left, double right);
std::vector<double> operator+(const std::vector<double>& left, const std::vector<double>& right);
std::vector<double> operator-(const std::vector<double>& left, const std::vector<double>& right);
double operator*(const std::vector<double>& left, const std::vector<double>& right);
bool equal(const std::vector<double>& left, const std::vector<double>& right, double eps);
bool operator==(const std::vector<double>& left, const std::vector<double>& right);
double norm2(const std::vector<double>& vec);
std::vector<double> normalize(const std::vector<double>& vec);
std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec);

namespace _random {
    unsigned getUnsigned(unsigned min, unsigned max);
    double getDouble(double min, double max);
    DenseMatrix getDenseMatrix(unsigned rows, unsigned cols, double density = 1, double min_el = 0, double max_el = 1, bool integer = false);
    CSRMatrix getCSRMatrix(unsigned rows, unsigned cols, double density, double min_el = 0, double max_el = 1);
    std::vector<double> getVector(unsigned rows, double min_el = 0, double max_el = 1);
    CSRMatrix getTestMatrix(unsigned rows);
}

double calcMaxEigenvalue(const CSRMatrix& A, double precision);

namespace JacobiTools {
    std::vector<double> multiply(const CSRMatrix& csr, const std::vector<double>& v);
    CSRMatrix inverseDiagonal(const CSRMatrix& csr);
}

namespace GaussSeidelTools {
    std::vector<double> inverseDiagonal(const CSRMatrix& csr);
}

namespace AccGaussSeidelTools {
    std::vector<double> iterate(const CSRMatrix& A, const std::vector<double>& d, const std::vector<double>& b, const std::vector<double>& x);
}

namespace FPIATools {
    // Не уверен, что правильно юзаю констэкспр, извини заранее за кринж, код в hpp и глобальные переменные
    constexpr unsigned nroots = 128;

    constexpr std::array<size_t, nroots> calcPermutations() {
        std::array<size_t, nroots> permutations = { 0 };

        // Логарифм не констэкспр, поэтому считаю так
        unsigned lg = 0, tmp = 1;
        while (tmp < nroots) {
            tmp *= 2;
            lg++;
        }

        size_t step = nroots;
        for (size_t i = 0; i < lg; i++) {
            for (size_t j = 0; j < nroots; j += step) {
                permutations[j + step / 2] = nroots / (step / 2) - 1 - permutations[j];
            }
            step /= 2;
        }
        return permutations;
    }

    constexpr std::array<double, nroots> calcChebyshevRoots() {
        std::array<double, nroots> roots = { 0 };

        // Синус тоже не констэкспр, поэтому тейлор (n большое)
        // Если бы тригонометрия была констэкспр, то можно
        // было бы обойтись одним синусом (а не четырьмя)!

        // Отличие в 12 знаке
        double s_p_2n = M_PI / 2 / nroots - (M_PI / 2 / nroots) * (M_PI / 2 / nroots) * (M_PI / 2 / nroots) / 6;
        // Отличие в 10 знаке
        double c_p_2n = 1 - (M_PI / 2 / nroots) * (M_PI / 2 / nroots) / 2;

        double s_p_n = 2 * s_p_2n * c_p_2n;
        double c_p_n = c_p_2n * c_p_2n - s_p_2n * s_p_2n;
        
        double s = s_p_2n;
        roots[0] = c_p_2n;
        for (size_t i = 1; i < nroots / 2; i++) {
            roots[i] = roots[i - 1] * c_p_n - s * s_p_n;
            s = s * c_p_n + roots[i - 1] * s_p_n;
        }
        for (size_t i = nroots / 2; i < nroots; i++) {
            roots[i] = -roots[nroots - i - 1];
        }

        return roots;
    }

    constexpr std::array<size_t, nroots> permutations = calcPermutations();
    constexpr std::array<double, nroots> roots = calcChebyshevRoots();
    
}
