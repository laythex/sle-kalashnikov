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

    std::ofstream fileFPI("../tests/accelerationsComparison/FPI.txt", std::ios::out);
    std::ofstream fileJacobi("../tests/accelerationComparison/Jacobi.txt", std::ios::out);
    std::ofstream fileGS("../tests/accelerationComparison/GS.txt", std::ios::out);
    std::ofstream fileFPIA("../tests/accelerationComparison/FPIA.txt", std::ios::out);

    unsigned m = 128;
    std::vector<double> roots = FPIAcceleratedTools::calcChebyshevRoots(m);
    std::vector<size_t> permutations = FPIAcceleratedTools::calcPermutations(m); 

    for (unsigned n = 0; n <= max_rows; n += rows_step) {
        A = _random::getDiagonallyDominantCSRMatrix(n, density, 0, 1);
        x_real = _random::getVector(n);
        b = A * x_real;
        x0 = _random::getVector(n);

        double lambda = FPIAcceleratedTools::calcMaxEigenvalue(A, 1e-3);

        auto startFPI = std::chrono::high_resolution_clock::now(); 
        x_calc = solvers::fixedPointIteration(A, b, 1, x0, 1e-12);
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
        x_calc = solvers::FPIAccelerated(A, b, lambda, roots, permutations, x0, 1e-12);
        auto endFPIA = std::chrono::high_resolution_clock::now();
        auto nsecFPIA = endFPIA - startFPIA;
        fileFPIA << n << '\t' << nsecFPIA.count() << std::endl;

        std::cout << n << std::endl;
    }

    fileFPI.close();
    fileJacobi.close();
    fileGS.close();

    return 0;
}