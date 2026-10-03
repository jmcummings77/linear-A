/* Bounded libFuzzer target for the public C API and shared native solver core. */
#include "../../ports/c/matrix.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void check(matrix_status status) { assert(status == M_OK); }
static void same(const matrix *a, const matrix *b) {
    assert(a->rows == b->rows && a->cols == b->cols);
    for(size_t i=0;i<a->rows*a->cols;i++) assert(a->values[i] == b->values[i]);
}
int LLVMFuzzerTestOneInput(const uint8_t *data, size_t length) {
    if(length < 4) return 0;
    const size_t m=data[0]%6,k=data[1]%6,n=data[2]%6;
    double av[25],bv[25];
    for(size_t i=0;i<25;i++) {
        av[i]=((int)data[(4+2*i)%length]-128)/16.0;
        bv[i]=((int)data[(5+2*i)%length]-128)/16.0;
    }
    matrix a={0},b={0},out={0},t={0},roundtrip={0},snapshot={0};
    check(m_from_array(m,k,av,&a));check(m_from_array(k,n,bv,&b));check(m_copy(&a,&snapshot));
    check(m_multiply(&a,&b,&out));
    assert(out.rows==m && out.cols==n);
    for(size_t i=0;i<m;i++) for(size_t j=0;j<n;j++) {
        double exact=0;
        for(size_t p=0;p<k;p++) exact+=av[i*k+p]*bv[p*n+j];
        assert(out.values[i*n+j]==exact); /* bounded dyadics are exact in binary64 */
    }
    same(&a,&snapshot);
    check(m_transpose(&a,&t));check(m_transpose(&t,&roundtrip));same(&a,&roundtrip);
    m_free(&t);m_free(&roundtrip);m_free(&out);
    /* Exercise invalid shape and index boundaries without invalid pointers. */
    double value=91;
    assert(m_get(&a,m,k,&value)==M_INDEX);assert(value==91);
    if(m!=k || k!=n) assert(m_add(&a,&b,&out)==M_SHAPE);
    double bad=NAN;
    assert(m_from_array(1,1,&bad,&out)==M_NONFINITE);
    assert(out.values==NULL);
    m_free(&a);m_free(&b);m_free(&snapshot);

    /* Small, well-conditioned SPD solves and factor reuse, including empty. */
    const size_t size=data[3]%6;
    double spd[25]={0},rhs[10]={0};
    for(size_t i=0;i<size;i++) for(size_t j=0;j<size;j++) {
        for(size_t p=0;p<size;p++) spd[i*size+j]+=av[p*size+i]*av[p*size+j];
        if(i==j)spd[i*size+j]+=1;
        rhs[i*2]+=spd[i*size+j];rhs[i*2+1]+=spd[i*size+j]*(double)(j+1);
    }
    check(m_from_array(size,size,spd,&a));check(m_from_array(size,2,rhs,&b));
    for(int algorithm=M_LU;algorithm<=M_QR;algorithm++) {
        matrix_factor *factor=NULL;
        check(m_factorize(&a,(matrix_factor_algorithm)algorithm,&factor));
        for(int repeat=0;repeat<2;repeat++) {
            check(m_factor_solve(factor,&b,&out));
            for(size_t i=0;i<size;i++) {
                assert(fabs(out.values[i*2]-1)<1e-9);
                assert(fabs(out.values[i*2+1]-(double)(i+1))<1e-9);
            }
            m_free(&out);
        }
        m_factor_free(&factor);assert(factor==NULL);
    }
    m_free(&a);m_free(&b);
    return 0;
}
