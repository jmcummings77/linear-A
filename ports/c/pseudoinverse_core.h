#ifndef LINEAR_A_PSEUDOINVERSE_CORE_H
#define LINEAR_A_PSEUDOINVERSE_CORE_H
#include <math.h>
#include <stddef.h>
/* Apply already-computed economy factors. A null RHS requests the inverse.
 * Returns 1 for invalid cutoff, 3 for range failure. Output is scratch storage. */
static inline size_t la_spectral_rank(size_t k,const double *s,double cutoff) {
    size_t rank=0;
    while(rank<k && s[rank]>0 && (cutoff==0 || s[rank]/s[0]>cutoff)) ++rank;
    return rank;
}
static inline double la_inverse_product(double a,double b,double c) {
    int ae,be,ce;
    double f=frexp(a,&ae)/frexp(b,&be)*frexp(c,&ce);
    return scalbn(f,ae-be+ce);
}
static inline int la_apply_inverse(size_t m,size_t n,size_t k,const double *u,
    const double *s,const double *vt,double cutoff,const double *rhs,size_t cols,
    int solving,double *out) {
    if(!isfinite(cutoff)||cutoff<0||cutoff>1)return 1;
    size_t rank=la_spectral_rank(k,s,cutoff);
    for(size_t j=0;j<cols;++j) {
        double scale=solving?0:1;
        if(solving)for(size_t i=0;i<m;++i)scale=fmax(scale,fabs(rhs[i*cols+j]));
        if(solving&&rank&&scale)for(size_t i=0;i<m;++i)
            if(rhs[i*cols+j]!=0 && rhs[i*cols+j]/scale==0)return 3;
        for(size_t p=0;p<rank;++p) {
            double projection=0;
            if(solving){if(scale)for(size_t i=0;i<m;++i)projection+=u[i*k+p]*(rhs[i*cols+j]/scale);}
            else projection=u[j*k+p];
            double coefficient=solving?la_inverse_product(projection,s[p],scale):0;
            if(!isfinite(coefficient))return 3;
            for(size_t i=0;i<n;++i){
                double term=solving?vt[p*n+i]*coefficient:la_inverse_product(projection,s[p],vt[p*n+i]);
                out[i*cols+j]+=term;
                if(!isfinite(out[i*cols+j]))return 3;
            }
        }
    }
    return 0;
}
#endif
