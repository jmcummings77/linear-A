#include "matrix.hpp"
#include <cmath>
#include <iostream>
int main() {
    linear_a::Matrix a(2,2,{4,1,2,3}), b(2,1,{6,8});
    auto x=a.solve(b), reconstructed=a*x;
    for (size_t i=0;i<2;i++)
        if (!std::isfinite(x.values()[i]) || std::abs(x.values()[i]-(i+1))>1e-12
            || !std::isfinite(reconstructed.values()[i]) || std::abs(reconstructed.values()[i]-b.values()[i])>1e-12) return 1;
    std::cout << "solution: 1, 2\n";
}
