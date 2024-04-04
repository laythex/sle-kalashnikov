#include "CSRMatrix.hpp"
#include "solvers.hpp"
#include "tools.hpp"

#include <iostream>
#include <fstream>
#include <chrono>

// Трехдиагональная матрица с диагнальным преобладанием
CSRMatrix getTestMatrix(unsigned rows) {
    std::vector<double> data(rows * rows);
    double a = _random::getDouble(0, 1);
    double b = _random::getDouble(0, a / 2);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < rows; j++) {
            if (i == j) {
                data[i * rows + j] = a;
            } else if (i == j + 1 || j == i + 1) {
                data[i * rows + j] = b;
            }
        }
    }

    return CSRMatrix(data, rows);
}

int main() {
    unsigned rows_step = 10, max_rows = 1000;

    CSRMatrix A;
    std::vector<double> b, x_real, x_calc, x0;

    std::ofstream fileFPI("../tests/accelerationsComparison/FPI.txt", std::ios::out);
    std::ofstream fileJacobi("../tests/accelerationComparison/Jacobi.txt", std::ios::out);
    std::ofstream fileGS("../tests/accelerationComparison/GS.txt", std::ios::out);
    std::ofstream fileFPIA("../tests/accelerationComparison/FPIA.txt", std::ios::out);

    for (unsigned n = 0; n <= max_rows; n += rows_step) {
        A = _random::getTestMatrix(n);
        x_real = _random::getVector(n);
        b = A * x_real;
        x0 = _random::getVector(n);

        double lambda = calcMaxEigenvalue(A, 1e-3);

        auto startFPI = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::fixedPointIteration(A, b, 1 / lambda, x0, 1e-12);
        auto endFPI = std::chrono::high_resolution_clock::now();
        auto nsecFPI = endFPI - startFPI;
        fileFPI << n << '\t' << nsecFPI.count() << std::endl;
    
        auto startJacobi = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::Jacobi(A, b, x0, 1e-12);
        auto endJacobi = std::chrono::high_resolution_clock::now();
        auto nsecJacobi = endJacobi - startJacobi;
        fileJacobi << n << '\t' << nsecJacobi.count() << std::endl;

        auto startGS = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::GaussSeidel(A, b, x0, 1e-12);
        auto endGS = std::chrono::high_resolution_clock::now();
        auto nsecGS = endGS - startGS;
        fileGS << n << '\t' << nsecGS.count() << std::endl;

        auto startFPIA = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::FPIAccelerated(A, b, 0, lambda, x0, 1e-12);
        auto endFPIA = std::chrono::high_resolution_clock::now();
        auto nsecFPIA = endFPIA - startFPIA;
        fileFPIA << n << '\t' << nsecFPIA.count() << std::endl;

        std::cout << n << std::endl;
    }

    fileFPI.close();
    fileJacobi.close();
    fileGS.close();
    fileFPIA.close();

    return 0;
}