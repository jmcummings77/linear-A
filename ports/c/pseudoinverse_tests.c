#include "matrix.h"
#include <assert.h>
#include <math.h>
int main(void){
    matrix a={0},b={0},p={0},x={0};double av[]={4,0,0,1},bv[]={4,8,1,2};
    assert(m_from_array(2,2,av,&a)==M_OK);assert(m_from_array(2,2,bv,&b)==M_OK);
    assert(m_pseudoinverse(&a,.25,&p)==M_OK);assert(p.values[0]==.25&&p.values[3]==0);
    assert(m_solve_minimum_norm(&a,&b,.25,&x)==M_OK);assert(x.values[0]==1&&x.values[1]==2&&x.values[2]==0&&x.values[3]==0);
    assert(m_pseudoinverse(&a,0,&a)==M_ARGUMENT);assert(m_pseudoinverse(&a,0,&p)==M_ARGUMENT);
    assert(m_solve_minimum_norm(&a,&b,0,&b)==M_ARGUMENT);
    matrix_spectral_diagnostics d={99,99,99};
    assert(m_spectral_diagnostics(&a,NAN,&d)==M_ARGUMENT);assert(d.rank==99);
    assert(m_spectral_diagnostics(&a,.25,&d)==M_OK);assert(d.rank==1&&d.reciprocal_condition==.25&&d.retained_reciprocal_condition==1);
    m_free(&x);m_free(&p);m_free(&a);m_free(&b);
    double tiny[]={1e-310};assert(m_from_array(1,1,tiny,&a)==M_OK);assert(m_from_array(1,1,tiny,&b)==M_OK);
    assert(m_pseudoinverse(&a,0,&p)==M_SOLVER_RANGE);assert(!p.values&&!p.rows&&!p.cols);
    assert(m_solve_minimum_norm(&a,&b,0,&x)==M_OK);assert(fabs(x.values[0]-1)<1e-14);
    m_free(&a);m_free(&b);m_free(&x);
    return 0;
}
