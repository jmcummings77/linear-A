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
    auto ridge=a.solve_ridge(b,1);
    if(std::abs(ridge.values()[0]-140.0/131)>1e-12||std::abs(ridge.values()[1]-230.0/131)>1e-12)return 1;
    linear_a::CSRMatrix sparse(2,2,{0,2,4},{0,1,0,1},{4,1,1,3});
    linear_a::ILU0 ilu(sparse);
    if(std::abs(ilu.apply({6,7})[0]-1)>1e-12||std::abs(ilu.apply({11,13})[0]-20.0/11)>1e-12||!sparse.gmres({11,13},20,1e-10,0,1000,false,false,&ilu).converged)return 1;
    auto gm=sparse.gmres({6,7},2,1e-10,0,1000,true,true);
    if(!gm.converged||std::abs(gm.x[0]-1)>1e-12||std::abs(gm.x[1]-2)>1e-12)return 1;
    auto product=sparse.matvec({1,2});auto cg=sparse.conjugate_gradient({6,7},1e-10,0,1000,true,true);
    if(product[0]!=6||product[1]!=7||!cg.converged||std::abs(cg.x[0]-1)>1e-12||std::abs(cg.x[1]-2)>1e-12)return 1;
    auto order=sparse.reverse_cuthill_mckee();auto reordered=sparse.permute_symmetric(order);auto y=linear_a::CSRMatrix::permute_vector(order,{1,2});
    if(linear_a::CSRMatrix::permute_vector(order,reordered.matvec(y),true)!=std::vector<double>({6,7}))return 1;
    if(sparse.approximate_minimum_degree()!=std::vector<std::size_t>{0,1})return 1;
    linear_a::SparseCholeskySymbolic plan(sparse);auto chol=plan.factorize(sparse);
    if(std::abs(chol.solve({6,7})[0]-1)>1e-12||std::abs(chol.solve({11,13})[0]-20.0/11)>1e-12||chol.lower().nnz()!=3)return 1;
    std::cout << "solution: 1, 2\n";
}
