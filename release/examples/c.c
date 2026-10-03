#include "matrix.h"
#include <math.h>
#include <stdio.h>
int main(void) {
    matrix a={0}, b={0}, x={0}, reconstructed={0};
    double av[]={4,1,2,3}, bv[]={6,8};
    int failed=m_from_array(2,2,av,&a) || m_from_array(2,1,bv,&b)
        || m_solve(&a,&b,&x) || m_multiply(&a,&x,&reconstructed);
    if (!failed) for (size_t i=0;i<2;i++)
        failed |= !isfinite(x.values[i]) || fabs(x.values[i]-(i+1))>1e-12
            || !isfinite(reconstructed.values[i]) || fabs(reconstructed.values[i]-bv[i])>1e-12;
    matrix u={0},s={0},vt={0};
    if (!failed) failed=m_svd(&a,1e-12,100,&u,&s,&vt);
    if (!failed) failed=s.rows!=2 || s.values[1]<=0;
    m_free(&u);m_free(&s);m_free(&vt);
    m_free(&reconstructed);m_free(&x);m_free(&b);m_free(&a);
    if (!failed) puts("solution: 1, 2");
    return failed;
}
