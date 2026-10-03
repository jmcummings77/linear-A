#ifndef LINEAR_A_CPP_SPARSE_HPP
#define LINEAR_A_CPP_SPARSE_HPP
#include "matrix.hpp"
#include "../c/sparse_core.h"
#include "../c/gmres_core.h"
namespace linear_a {
struct CGResult {std::vector<double> x,residuals;std::vector<std::vector<double>> iterates;std::size_t iterations;std::string reason;bool converged;};
struct GMRESResult {std::vector<double> x,residuals,estimated_residuals;std::vector<std::vector<double>> iterates;std::vector<std::size_t> restarts;std::size_t iterations;std::string reason;bool converged;};
class CSRMatrix;
class ILU0 {
 friend class CSRMatrix;
 la_ilu0 factor_{};
public:
 explicit ILU0(const CSRMatrix& a);
 ~ILU0(){la_ilu0_free(&factor_);}
 ILU0(const ILU0&)=delete;ILU0& operator=(const ILU0&)=delete;
 std::size_t size()const{return factor_.factors.rows;}std::size_t nnz()const{return factor_.factors.nnz;}
 std::vector<double> apply(const std::vector<double>& b)const{std::vector<double> x(size());int c=la_ilu0_apply(&factor_,b.data(),b.size(),x.data());if(c)throw std::runtime_error("invalid or nonfinite ILU0 solve");return x;}
};
class CSRMatrix {
 friend class ILU0;
    std::size_t rows_,cols_;std::vector<std::size_t> rp_,ci_;std::vector<double> v_;
    la_csr view() const {return {rows_,cols_,v_.size(),const_cast<std::size_t*>(rp_.data()),const_cast<std::size_t*>(ci_.data()),const_cast<double*>(v_.data())};}
public:
    CSRMatrix(std::size_t rows,std::size_t cols,std::vector<std::size_t> offsets,std::vector<std::size_t> indices,std::vector<double> values):rows_(rows),cols_(cols),rp_(std::move(offsets)),ci_(std::move(indices)),v_(std::move(values)) {
        if(rows==SIZE_MAX||rp_.size()!=rows+1||ci_.size()!=v_.size())throw std::invalid_argument("invalid CSR arrays");
        auto a=view();if(la_csr_validate(&a))throw std::invalid_argument("invalid canonical CSR");
    }
    static CSRMatrix from_dense(const Matrix& a){
        std::vector<std::size_t> rp{0},ci;std::vector<double> v;
        for(std::size_t i=0;i<a.rows();i++){for(std::size_t j=0;j<a.cols();j++)if(a.get(i,j)!=0){ci.push_back(j);v.push_back(a.get(i,j));}rp.push_back(v.size());}
        return {a.rows(),a.cols(),std::move(rp),std::move(ci),std::move(v)};
    }
    std::size_t rows()const{return rows_;}std::size_t cols()const{return cols_;}std::size_t nnz()const{return v_.size();}
    const std::vector<std::size_t>& row_offsets()const{return rp_;}const std::vector<std::size_t>& column_indices()const{return ci_;}const std::vector<double>& values()const{return v_;}
    std::vector<double> matvec(const std::vector<double>& x)const{
        if(x.size()!=cols_)throw std::invalid_argument("incompatible vector");
        auto a=view();std::vector<double> out(rows_);
        if(la_csr_mv(&a,x.data(),out.data()))throw std::overflow_error("sparse multiplication outside float64 range");return out;
    }
    CGResult conjugate_gradient(const std::vector<double>& b,double rtol=1e-10,double atol=0,std::size_t limit=1000,bool jacobi=false,bool capture=false)const{
        auto a=view();la_cg_result r={};int code=la_csr_cg(&a,b.data(),b.size(),rtol,atol,limit,jacobi,capture,&r);
        if(code==2)throw std::bad_alloc();if(code)throw std::invalid_argument("invalid CG input, options, symmetry or Jacobi diagonal");
        const char *reasons[]={"converged","iteration_limit","breakdown","nonfinite"};
        try{
            CGResult result;result.x.assign(r.x,r.x+r.size);result.residuals.assign(r.residuals,r.residuals+r.iterations+1);result.iterations=r.iterations;result.reason=reasons[r.reason];result.converged=r.reason==0;
            if(capture)for(std::size_t i=0;i<=r.iterations;i++)result.iterates.emplace_back(r.iterates+i*r.size,r.iterates+(i+1)*r.size);
            la_cg_free(&r);return result;
        }catch(...){la_cg_free(&r);throw;}
    }
    GMRESResult gmres(const std::vector<double>& b,std::size_t restart=30,double rtol=1e-10,double atol=0,std::size_t limit=1000,bool jacobi=false,bool capture=false,const ILU0 *preconditioner=nullptr)const{
        auto a=view();la_gmres_result r={};int code=la_csr_gmres_preconditioned(&a,b.data(),b.size(),restart,rtol,atol,limit,jacobi,capture,preconditioner?&preconditioner->factor_:nullptr,&r);
        if(code==2)throw std::bad_alloc();if(code)throw std::invalid_argument("invalid GMRES input, options or Jacobi diagonal");
        const char *reasons[]={"converged","iteration_limit","breakdown","nonfinite","stagnation"};
        try{GMRESResult result;result.x.assign(r.x,r.x+r.size);result.residuals.assign(r.residuals,r.residuals+r.iterations+1);result.estimated_residuals.assign(r.estimated_residuals,r.estimated_residuals+r.iterations+1);result.restarts.assign(r.restarts,r.restarts+r.restart_count);result.iterations=r.iterations;result.reason=reasons[r.reason];result.converged=r.reason==0;
            if(capture)for(std::size_t i=0;i<=r.iterations;i++)result.iterates.emplace_back(r.iterates+i*r.size,r.iterates+(i+1)*r.size);
            la_gmres_free(&r);return result;
        }catch(...){la_gmres_free(&r);throw;}
    }

};
inline ILU0::ILU0(const CSRMatrix& a){auto view=a.view();int c=la_ilu0_create(&view,&factor_);if(c==2)throw std::bad_alloc();if(c)throw std::runtime_error("invalid ILU0 matrix, zero pivot or nonfinite factor");}
}
#endif
