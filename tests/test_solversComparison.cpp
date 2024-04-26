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

    std::ofstream fileFPI("../tests/solversComparison/FPI.txt", std::ios::out);
    std::ofstream fileFPIA("../tests/solversComparison/FPIA.txt", std::ios::out);
    std::ofstream fileJacobi("../tests/solversComparison/Jacobi.txt", std::ios::out);
    std::ofstream fileGS("../tests/solversComparison/GS.txt", std::ios::out);
    std::ofstream fileSGS("../tests/solversComparison/SGS.txt", std::ios::out);
    std::ofstream fileAGS("../tests/solversComparison/AGS.txt", std::ios::out);
    std::ofstream fileGradDes("../tests/solversComparison/GradDes.txt", std::ios::out);
    std::ofstream fileConjGrad("../tests/solversComparison/ConjGrad.txt", std::ios::out);
    std::ofstream fileGMRES("../tests/solversComparison/GMRES.txt", std::ios::out);

    for (unsigned n = 0; n <= max_rows; n += rows_step) {
        A = _random::getTestMatrix(n);
        x_real = _random::getVector(n);
        b = A * x_real;
        x0 = _random::getVector(n);

        double lambda = calcMaxEigenvalue(A, 1e-3);

        auto start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::fixedPointIteration(A, b, 1 / lambda, x0, 1e-12);
        auto end = std::chrono::high_resolution_clock::now();
        auto nsec = end - start;
        fileFPI << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::FPIAccelerated(A, b, 0, lambda, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileFPIA << n << '\t' << nsec.count() << std::endl;
    
        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::Jacobi(A, b, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileJacobi << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::GaussSeidel(A, b, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileGS << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::SymGaussSeidel(A, b, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileSGS << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::AccGaussSeidel(A, b, 0.9, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileAGS << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::GradientDescent(A, b, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileGradDes << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::GradientDescent(A, b, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileConjGrad << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::GradientDescent(A, b, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileConjGrad << n << '\t' << nsec.count() << std::endl;

        start = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::GMRES(A, b, x0, 1e-12);
        end = std::chrono::high_resolution_clock::now();
        nsec = end - start;
        fileGMRES << n << '\t' << nsec.count() << std::endl;

        std::cout << n << std::endl;
    }

    fileFPI.close();
    fileFPIA.close();
    fileJacobi.close();
    fileGS.close();
    fileSGS.close();
    fileAGS.close();
    fileGradDes.close();
    fileConjGrad.close();
    fileGMRES.close();

    return 0;
}