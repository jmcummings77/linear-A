#include "matrix.h"
#include <assert.h>
#include <math.h>
int main(void){
 size_t rp[]={0,2,4},ci[]={0,1,0,1};double v[]={4,1,1,3},b[]={6,7};sparse_matrix a={0};matrix_ic0 f={0};matrix_cg_result r={0};
 assert(m_csr_create(2,2,4,rp,ci,v,&a)==M_OK);assert(m_ic0_create(&a,&f)==M_OK);assert(f.lower.nnz==3);
 assert(m_csr_cg_preconditioned(&a,b,2,1e-10,0,100,false,true,&f,&r)==M_OK);assert(r.reason==0&&r.iterations==1);assert(fabs(r.x[0]-1)<1e-12&&fabs(r.x[1]-2)<1e-12);m_cg_free(&r);
 assert(m_csr_cg_preconditioned(&a,b,2,1e-10,0,100,true,false,&f,&r)==M_ARGUMENT);assert(!r.x);
 assert(m_ic0_create(&a,&f)==M_ARGUMENT);m_ic0_free(&f);m_ic0_free(&f);m_csr_free(&a);
 size_t sp[]={0,3,6,9,12},sc[]={0,1,3,0,1,2,1,2,3,0,2,3};double sv[]={1.5,1,1,1,1.5,1,1,1.5,-1,1,-1,1.5};assert(m_csr_create(4,4,12,sp,sc,sv,&a)==M_OK);
 assert(m_ic0_create(&a,&f)==M_NOT_POSITIVE_DEFINITE);assert(!f.lower.values);matrix_cholesky_symbolic s={0};matrix_cholesky full={0};assert(m_cholesky_analyze(&a,&s)==M_OK);assert(m_cholesky_factorize(&s,&a,&full)==M_OK);m_cholesky_free(&full);m_cholesky_symbolic_free(&s);m_csr_free(&a);
 return 0;
}
