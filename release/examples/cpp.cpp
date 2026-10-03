#include "matrix.hpp"
#include <cmath>
#include <iostream>
int main() {
    linear_a::Matrix a(2,2,{4,1,2,3}), b(2,1,{6,8});
    auto x=a.solve(b), reconstructed=a*x;
    for (size_t i=0;i<2;i++)
        if (!std::isfinite(x.values()[i]) || std::abs(x.values()[i]-(i+1))>1e-12
            || !std::isfinite(reconstructed.values()[i]) || std::abs(reconstructed.values()[i]-b.values()[i])>1e-12) return 1;
    auto decomposition=a.svd();
    if(decomposition.values.size()!=2 || decomposition.values[1]<=0) return 1;
    auto inverse=a.pseudoinverse(),minimum=a.solve_minimum_norm(b);
    if(std::abs(inverse.values()[0]-.3)>1e-12 || std::abs(minimum.values()[1]-2)>1e-12 || a.spectral_diagnostics().rank!=2)return 1;
    std::cout << "solution: 1, 2\n";
}
