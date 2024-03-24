#pragma once

#include <vector>
#include <iostream>

class TridiagMatrix {
    public:
        void initialize(const std::vector<double>& a_short_data, const std::vector<double>& b_data, const std::vector<double>& c_short_data);

        std::vector<double> solve(const std::vector<double>& d) const;
    
    protected:
        unsigned n;
        std::vector<double> a, b, c;
};
