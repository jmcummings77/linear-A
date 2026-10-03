#include "matrix.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
int main(void){
    matrix a={0},u={0},s={0},vt={0};double data[]={1,2,3,4,5,7,6,8,9};
    assert(m_from_array(3,3,data,&a)==M_OK);
    assert(m_svd(&a,0,100,&u,&s,&vt)==M_ARGUMENT);
    assert(m_svd(&a,1e-12,100,&u,&u,&vt)==M_ARGUMENT);
    assert(m_svd(&a,1e-12,100,&a,&s,&vt)==M_ARGUMENT);
    assert(m_create(1,1,&u)==M_OK);
    assert(m_svd(&a,1e-12,100,&u,&s,&vt)==M_ARGUMENT);m_free(&u);
    for(int i=0;i<100;i++){
        assert(m_svd(&a,1e-12,1,&u,&s,&vt)==M_NO_CONVERGENCE);
        assert(!u.values&&!s.values&&!vt.values&&!u.rows&&!s.rows&&!vt.rows);
        assert(m_svd(&a,1e-12,100,&u,&s,&vt)==M_OK);
        assert(a.values[0]==1);m_free(&u);m_free(&s);m_free(&vt);
    }
    a.values[0]=NAN;assert(m_svd(&a,1e-12,100,&u,&s,&vt)==M_NONFINITE);m_free(&a);
    matrix bad={1,1,NULL};assert(m_svd(&bad,1e-12,100,&u,&s,&vt)==M_ARGUMENT);
    matrix empty={SIZE_MAX,0,NULL};assert(m_svd(&empty,1e-12,100,&u,&s,&vt)==M_OK);
    m_free(&u);m_free(&s);m_free(&vt);
    return 0;
}
