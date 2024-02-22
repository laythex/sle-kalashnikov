#include "CSRMatrix.hpp"
#include "DenseMatrix.hpp"
#include "tools.hpp"

#include <iostream>
#include <fstream>
#include <chrono>
#include <fstream>

int main() {
    std::ofstream file("../out.txt", std::ios::out);

    unsigned cols_step = 20, max_cols = 1000;
    double density_step = 0.005;
    std::vector<double> vec, res_dense, res_sparse;
    CSRMatrix sparse;
    DenseMatrix dense;

    file << cols_step << '\t' << max_cols << '\t' << density_step << std::endl;

    for (unsigned cols = cols_step; cols <= max_cols; cols += cols_step) {
        vec = randomVector(cols);
        dense = randomDenseMatrix(cols, cols);

        std::cout << cols << std::endl;

        for (double density = density_step; density <= 1; density += density_step) {
            sparse = randomCSRMatrix(cols, cols, density);

            auto start_dense = std::chrono::high_resolution_clock::now(); 
            res_dense = dense * vec;
            auto end_dense = std::chrono::high_resolution_clock::now();
            auto nsec_dense = end_dense - start_dense;

            auto start_sparse = std::chrono::high_resolution_clock::now();
            res_sparse = sparse * vec;
            auto end_sparse = std::chrono::high_resolution_clock::now();
            auto nsec_sparse = end_sparse - start_sparse;
            
            double factor = static_cast<double>(nsec_dense.count() / nsec_sparse.count());
            file << factor << '\t';

            // Чтобы оптимизация не съела умножения
            if (res_dense.size() != res_sparse.size()) {
                std::cout << "Failure!";
                return 0;
            }
        }

        file << std::endl;
    }

    return 0;
}