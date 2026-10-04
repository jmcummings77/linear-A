#include "matrix.h"
#include <assert.h>
#include <math.h>
static int reject_frame(void *context,size_t level,size_t width,int phase,const double*x,const double*b,const double*r){(void)context;(void)level;(void)width;(void)phase;(void)x;(void)b;(void)r;return 2;}
int main(void){
 matrix_multigrid m={0};assert(m_multigrid_create(2,&m)==M_ARGUMENT);assert(!m.width);assert(m_multigrid_create(7,&m)==M_OK);assert(m_multigrid_create(7,&m)==M_ARGUMENT);
 sparse_matrix a={0};matrix b={0},x={0};assert(m_multigrid_matrix(&m,&a)==M_OK);assert(a.rows==49&&a.nnz==217);assert(m_create(49,1,&b)==M_OK);for(size_t i=0;i<49;i++)b.values[i]=sin((double)i+1);
 assert(m_multigrid_apply(&m,&b,&x)==M_OK);double first=x.values[0];m_free(&x);assert(m_multigrid_apply(&m,&b,&x)==M_OK);assert(x.values[0]==first);m_free(&x);
 double target[49];for(size_t i=0;i<49;i++)target[i]=99;assert(la_mg_apply_trace(&m,b.values,49,target,reject_frame,NULL)==2);for(size_t i=0;i<49;i++)assert(target[i]==99);
 matrix_cg_result r={0};assert(m_csr_cg_multigrid(&a,b.values,49,1e-10,0,100,false,true,&m,&r)==M_OK);assert(r.reason==0&&r.iterations<20);m_cg_free(&r);
 assert(m_csr_cg_multigrid(&a,b.values,49,1e-10,0,100,true,false,&m,&r)==M_ARGUMENT);assert(!r.x);
 b.values[0]=NAN;assert(m_multigrid_apply(&m,&b,&x)!=M_OK);assert(!x.values);
 for(size_t i=0;i<49;i++)b.values[i]=1.7e308;assert(m_multigrid_apply(&m,&b,&x)==M_SOLVER_RANGE);assert(!x.values);
 m_free(&b);m_csr_free(&a);return 0;
}
