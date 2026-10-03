#include "matrix.h"
#include <assert.h>
#include <math.h>
int main(void){
 size_t rp[]={0,2,4},ci[]={0,1,0,1};double v[]={4,1,1,3},b[]={6,7};sparse_matrix a={0};
 assert(m_csr_create(2,2,4,rp,ci,v,&a)==M_OK);v[0]=99;assert(a.values[0]==4);
 matrix_cg_result r={0};assert(m_csr_cg(&a,b,2,1e-10,0,100,1,1,&r)==M_OK);
 assert(r.reason==0&&r.iterations==2&&fabs(r.x[0]-1)<1e-12&&fabs(r.x[1]-2)<1e-12);
 assert(r.iterates[0]==0&&r.iterates[1]==0);assert(r.residuals[r.iterations]<1e-10);
 assert(m_csr_cg(&a,b,2,1e-10,0,100,0,0,&r)==M_ARGUMENT);m_cg_free(&r);m_cg_free(&r);
 matrix x={0},out={0};assert(m_create(2,1,&x)==M_OK);x.values[0]=1;x.values[1]=2;
 assert(m_csr_matvec(&a,&x,&x)==M_ARGUMENT);assert(m_csr_matvec(&a,&x,&out)==M_OK);assert(out.values[0]==6&&out.values[1]==7);m_free(&x);m_free(&out);
 assert(m_csr_create(2,2,4,rp,ci,v,&a)==M_ARGUMENT);m_csr_free(&a);m_csr_free(&a);
 ci[1]=0;assert(m_csr_create(2,2,4,rp,ci,v,&a)==M_ARGUMENT);assert(!a.offsets);
 size_t empty[]={0};assert(m_csr_create(0,0,0,empty,NULL,NULL,&a)==M_OK);assert(m_csr_cg(&a,NULL,0,1e-10,0,0,0,1,&r)==M_OK);assert(r.reason==0&&r.iterations==0);m_cg_free(&r);m_csr_free(&a);
 return 0;
}
