#include "CSR.hpp"
#include "DenseMatrix.hpp"
#include <iostream>

int main() {

    std::vector<double> data1 = {1, 2, 0, 0, 4, 0, 0, 2, 6};
    std::vector<double> data2 = {1, 2, 0, 3, 0, 0, 4, 0, 0, 1, 0, 11};

    CSR csr1, csr2;
    csr1.initialize(data1, 3);
    csr2.initialize(data2, 4);

    DenseMatrix dm1, dm2;
    dm1.initialize(data1, 3);
    dm2.initialize(data2, 4);

    std::cout << csr1 << std::endl;
    std::cout << csr2 << std::endl;

    std::cout << dm1 << std::endl;
    std::cout << dm2 << std::endl;

    return 0;
}