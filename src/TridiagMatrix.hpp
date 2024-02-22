#pragma once

#include <vector>
#include <iostream>

class TridiagMatrix {
    public:
        void initialize(const std::vector<double>& a_short_data, const std::vector<double>& b_data, const std::vector<double>& c_short_data);
        
        // double operator()(int i, int j) const;
        // std::vector<double> operator*(const std::vector<double>& v) const;
        // TridiagMatrix operator+(const TridiagMatrix& other) const;
        // double operator*(const TridiagMatrix& other) const;
        // friend std::ostream& operator<<(std::ostream& os, const TridiagMatrix& tm);

        std::vector<double> solve(const std::vector<double>& d) const;
    
    protected:
        unsigned n;
        std::vector<double> a, b, c;
};
