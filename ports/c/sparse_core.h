#ifndef LINEAR_A_SPARSE_CORE_H
#define LINEAR_A_SPARSE_CORE_H
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
/* Shared C/C++ numerical kernel. Codes: 1 argument, 2 allocation,
 * 3 nonfinite arithmetic, 4 asymmetry, 5 invalid Jacobi diagonal. */
typedef struct {size_t rows,cols,nnz;size_t *offsets,*indices;double *values;} la_csr;
typedef struct {double *x,*residuals,*iterates;size_t size,iterations;int reason;} la_cg_result;
static inline int la_sparse_count(size_t n,size_t width){return n<=SIZE_MAX/width;}
static inline void la_csr_free(la_csr *a){if(a){free(a->offsets);free(a->indices);free(a->values);memset(a,0,sizeof(*a));}}
static inline void la_cg_free(la_cg_result *r){if(r){free(r->x);free(r->residuals);free(r->iterates);memset(r,0,sizeof(*r));}}
static inline int la_csr_validate(const la_csr *a){
    if(!a||a->rows==SIZE_MAX||!a->offsets||(!a->indices&&a->nnz)||(!a->values&&a->nnz)||a->offsets[0]!=0||a->offsets[a->rows]!=a->nnz)return 1;
    for(size_t i=0;i<a->rows;i++){
        if(a->offsets[i]>a->offsets[i+1]||a->offsets[i+1]>a->nnz)return 1;
        for(size_t p=a->offsets[i];p<a->offsets[i+1];p++)
            if(a->indices[p]>=a->cols||!isfinite(a->values[p])||(p>a->offsets[i]&&a->indices[p]<=a->indices[p-1]))return 1;
    }
    return 0;
}
static inline int la_csr_create(size_t rows,size_t cols,size_t nnz,const size_t *rp,const size_t *ci,const double *v,la_csr *out){
    if(!out||out->offsets||out->indices||out->values||out->rows||out->cols||out->nnz||rows==SIZE_MAX)return 1;
    la_csr a={rows,cols,nnz,(size_t*)rp,(size_t*)ci,(double*)v};
    if(!la_sparse_count(rows+1,sizeof(size_t))||!la_sparse_count(nnz,sizeof(size_t))||!la_sparse_count(nnz,sizeof(double)))return 2;
    int code=la_csr_validate(&a);if(code)return code;
    a.offsets=(size_t*)malloc((rows+1)*sizeof(size_t));a.indices=(size_t*)malloc((nnz?nnz:1)*sizeof(size_t));a.values=(double*)malloc((nnz?nnz:1)*sizeof(double));
    if(!a.offsets||!a.indices||!a.values){la_csr_free(&a);return 2;}
    memcpy(a.offsets,rp,(rows+1)*sizeof(size_t));if(nnz){memcpy(a.indices,ci,nnz*sizeof(size_t));memcpy(a.values,v,nnz*sizeof(double));}
    *out=a;return 0;
}
static inline size_t la_csr_find(const la_csr *a,size_t row,size_t col){
    size_t lo=a->offsets[row],hi=a->offsets[row+1];
    while(lo<hi){size_t mid=lo+(hi-lo)/2;if(a->indices[mid]<col)lo=mid+1;else hi=mid;}
    return lo;
}
static inline int la_csr_mv(const la_csr *a,const double *x,double *out){
    for(size_t j=0;j<a->cols;j++)if(!isfinite(x[j]))return 3;
    for(size_t i=0;i<a->rows;i++){
        double sum=0;for(size_t p=a->offsets[i];p<a->offsets[i+1];p++)sum+=a->values[p]*x[a->indices[p]];
        if(!isfinite(sum))return 3;out[i]=sum;
    }return 0;
}
static inline double la_sparse_norm(size_t n,const double *x){double s=0;for(size_t i=0;i<n;i++)s=hypot(s,x[i]);return s;}
static inline double la_sparse_dot(size_t n,const double *x,const double *y){double s=0;for(size_t i=0;i<n;i++)s+=x[i]*y[i];return s;}
static inline int la_csr_cg(const la_csr *a,const double *b,size_t count,double rtol,double atol,size_t limit,int jacobi,int capture,la_cg_result *out){
    if(!out||out->x||out->residuals||out->iterates||out->size||out->iterations||out->reason)return 1;
    int code=la_csr_validate(a);if(code)return code;
    size_t n=a->rows;
    if(n!=a->cols||count!=n||(!b&&n)||!isfinite(rtol)||rtol<0||rtol>=1||!isfinite(atol)||atol<0||limit>100000||(jacobi!=0&&jacobi!=1)||(capture!=0&&capture!=1))return 1;
    for(size_t i=0;i<n;i++)if(!isfinite(b[i]))return 1;
    if(!la_sparse_count(n,8*sizeof(double))||(capture&&n>SIZE_MAX/(limit+1)/sizeof(double)))return 2;
    double *work=(double*)calloc(n?n*8:1,sizeof(double));
    if(!work)return 2;
    double *r=work,*z=r+n,*p=z+n,*q=p+n,*candidate=q+n,*res=candidate+n,*diag=res+n,*ax=diag+n;
    for(size_t i=0;i<n;i++){
        for(size_t k=a->offsets[i];k<a->offsets[i+1];k++){
            size_t j=a->indices[k],s=la_csr_find(a,j,i);
            double other=s<a->offsets[j+1]&&a->indices[s]==i?a->values[s]:0;
            if(a->values[k]!=other){free(work);return 4;}
        }
        diag[i]=1;
        if(jacobi){size_t k=la_csr_find(a,i,i);if(k==a->offsets[i+1]||a->indices[k]!=i||a->values[k]<=0){free(work);return 5;}diag[i]=a->values[k];}
    }
    la_cg_result result={NULL,NULL,NULL,n,0,1};
    result.x=(double*)calloc(n?n:1,sizeof(double));result.residuals=(double*)calloc(limit+1,sizeof(double));
    if(capture)result.iterates=(double*)calloc(n?(limit+1)*n:1,sizeof(double));
    if(!result.x||!result.residuals||(capture&&!result.iterates)){la_cg_free(&result);free(work);return 2;}
    if(n)memcpy(r,b,n*sizeof(double));result.residuals[0]=la_sparse_norm(n,r);
    double threshold=fmax(atol,rtol*result.residuals[0]),rho=0;
    if(!isfinite(result.residuals[0]))result.reason=3;
    else if(result.residuals[0]<=threshold)result.reason=0;
    else {
        for(size_t i=0;i<n;i++)p[i]=z[i]=r[i]/diag[i];rho=la_sparse_dot(n,r,z);
        for(size_t step=0;step<limit;step++){
            if(la_csr_mv(a,p,q)){result.reason=3;break;}
            double curvature=la_sparse_dot(n,p,q);
            if(!isfinite(rho)||!isfinite(curvature)){result.reason=3;break;}
            if(rho<=0||curvature<=0){result.reason=2;break;}
            double alpha=rho/curvature;for(size_t i=0;i<n;i++)candidate[i]=result.x[i]+alpha*p[i];
            if(la_csr_mv(a,candidate,ax)){result.reason=3;break;}
            for(size_t i=0;i<n;i++)res[i]=b[i]-ax[i];double length=la_sparse_norm(n,res);
            if(!isfinite(length)){result.reason=3;break;}
            if(n){memcpy(result.x,candidate,n*sizeof(double));memcpy(r,res,n*sizeof(double));}
            result.iterations++;result.residuals[result.iterations]=length;
            if(capture&&n)memcpy(result.iterates+result.iterations*n,result.x,n*sizeof(double));
            if(length<=threshold){result.reason=0;break;}
            for(size_t i=0;i<n;i++)z[i]=r[i]/diag[i];double next=la_sparse_dot(n,r,z);
            if(!isfinite(next)){result.reason=3;break;}if(next<=0){result.reason=2;break;}
            double beta=next/rho;for(size_t i=0;i<n;i++)p[i]=z[i]+beta*p[i];rho=next;
        }
    }
    free(work);*out=result;return 0;
}
#endif
