#include "CSRMatrix.hpp"
#include "solvers.hpp"
#include "tools.hpp"

#include <iostream>
#include <fstream>
#include <chrono>

int main() {
    unsigned rows_step = 10, max_rows = 1000;
    double density = 0.01;

    CSRMatrix A;
    std::vector<double> b, x_real, x_calc, x0;

    std::ofstream fileFPI("../tests/solversComparison/FPI.txt", std::ios::out);
    std::ofstream fileJacobi("../tests/solversComparison/Jacobi.txt", std::ios::out);
    std::ofstream fileGS("../tests/solversComparison/GS.txt", std::ios::out);

    for (unsigned n = 0; n <= max_rows; n += rows_step) {
        A = _random::getDiagonallyDominantCSRMatrix(n, density, 0, 1);
        x_real = _random::getVector(n);
        b = A * x_real;
        x0 = _random::getVector(n);

        auto startFPI = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::fixedPointIteration(A, b, 1, x0, 1e-14);
        auto endFPI = std::chrono::high_resolution_clock::now();
        auto nsecFPI = endFPI - startFPI;
        fileFPI << n << '\t' << nsecFPI.count() << std::endl;
    
        auto startJacobi = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::Jacobi(A, b, x0, 1e-14);
        auto endJacobi = std::chrono::high_resolution_clock::now();
        auto nsecJacobi = endJacobi - startJacobi;
        fileJacobi << n << '\t' << nsecJacobi.count() << std::endl;

        auto startGS = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::GaussSeidel(A, b, x0, 1e-14);
        auto endGS = std::chrono::high_resolution_clock::now();
        auto nsecGS = endGS - startGS;
        fileGS << n << '\t' << nsecGS.count() << std::endl;

        std::cout << n << std::endl;
    }

    fileFPI.close();
    fileJacobi.close();
    fileGS.close();

    return 0;
}