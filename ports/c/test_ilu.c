#include "matrix.h"
#include <assert.h>
#include <math.h>
int main(void){
 size_t rp[]={0,2,4},ci[]={0,1,0,1};double values[]={4,1,1,3};sparse_matrix a={0};matrix_ilu0 f={0};matrix b={0},x={0};
 assert(m_csr_create(2,2,4,rp,ci,values,&a)==M_OK);assert(m_ilu0_create(&a,&f)==M_OK);assert(m_ilu0_create(&a,&f)==M_ARGUMENT);m_csr_free(&a);
 assert(m_create(2,1,&b)==M_OK);b.values[0]=6;b.values[1]=7;assert(m_ilu0_apply(&f,&b,&x)==M_OK);assert(fabs(x.values[0]-1)<1e-12&&fabs(x.values[1]-2)<1e-12);m_free(&x);
 b.values[0]=11;b.values[1]=13;assert(m_ilu0_apply(&f,&b,&x)==M_OK);assert(fabs(x.values[0]-20./11)<1e-12&&fabs(x.values[1]-41./11)<1e-12);m_free(&x);m_free(&b);m_ilu0_free(&f);m_ilu0_free(&f);
 double singular[]={1,1,1,1};assert(m_csr_create(2,2,4,rp,ci,singular,&a)==M_OK);assert(m_ilu0_create(&a,&f)==M_SINGULAR);assert(!f.factors.values&&!f.diagonal);m_csr_free(&a);
 double overflow[]={1e-308,1e308,1e308,1};assert(m_csr_create(2,2,4,rp,ci,overflow,&a)==M_OK);assert(m_ilu0_create(&a,&f)==M_SOLVER_RANGE);assert(!f.factors.values&&!f.diagonal);m_csr_free(&a);
 size_t empty[]={0};assert(m_csr_create(0,0,0,empty,NULL,NULL,&a)==M_OK);assert(m_ilu0_create(&a,&f)==M_OK);m_ilu0_free(&f);m_csr_free(&a);
 return 0;
}
