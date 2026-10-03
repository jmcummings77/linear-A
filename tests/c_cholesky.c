#include "../ports/c/matrix.h"
#include <assert.h>
#include <math.h>
int main(void){
 size_t rp[]={0,2,4},ci[]={0,1,0,1};double v[]={4,2,2,3};sparse_matrix a={0};matrix_cholesky_symbolic s={0};matrix_cholesky f={0},g={0};
 assert(m_csr_create(2,2,4,rp,ci,v,&a)==M_OK);assert(m_cholesky_analyze(&a,&s)==M_OK);assert(m_cholesky_factorize(&s,&a,&f)==M_OK);
 for(size_t p=0;p<4;p++)a.values[p]*=2;assert(m_cholesky_factorize(&s,&a,&g)==M_OK);
 a.values[0]=-1;matrix_cholesky bad={0};assert(m_cholesky_factorize(&s,&a,&bad)==M_NOT_POSITIVE_DEFINITE);assert(!bad.lower.offsets);
 a.values[0]=8;a.values[1]=3;assert(m_cholesky_factorize(&s,&a,&bad)==M_NOT_SYMMETRIC);assert(!bad.lower.offsets);
 m_csr_free(&a);m_cholesky_symbolic_free(&s);m_cholesky_symbolic_free(&s);
 double rhs[]={8,8},x[2];assert(la_cholesky_solve(&f,rhs,2,x)==0);assert(fabs(x[0]-1)<1e-12&&fabs(x[1]-2)<1e-12);
 double twice[]={16,16};assert(la_cholesky_solve(&g,twice,2,x)==0);assert(fabs(x[0]-1)<1e-12&&fabs(x[1]-2)<1e-12);
 assert(la_cholesky_solve(&f,rhs,1,x)==1);rhs[0]=NAN;assert(la_cholesky_solve(&f,rhs,2,x)==1);
 rhs[0]=rhs[1]=0;assert(la_cholesky_solve(&f,rhs,2,x)==0&&x[0]==0&&x[1]==0);
 m_cholesky_free(&f);m_cholesky_free(&g);m_cholesky_free(&f);return 0;
}
