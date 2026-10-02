#ifndef LINEAR_A_SOLVE_CORE_H
#define LINEAR_A_SOLVE_CORE_H
/* Shared C/C++ factor arithmetic. Factors own a scaled snapshot of the input. */
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
enum { LA_LU=1, LA_CHOLESKY=2, LA_QR=3 };
enum { LA_OK, LA_SHAPE, LA_MEMORY, LA_RANGE, LA_SINGULAR, LA_NOT_SPD, LA_RANK, LA_ARGUMENT };
typedef struct {
    size_t rows, cols;
    int algorithm;
    double scale, norm, *data, *tau;
    size_t *permutation;
} la_factor;
static void la_factor_destroy(la_factor *f) {
    if(f){free(f->data);free(f->tau);free(f->permutation);free(f);}
}
static int la_scaled(double value,double scale,double *out) {
    *out=value/scale;
    return isfinite(*out)&&(!value||*out!=0)?LA_OK:LA_RANGE;
}
static int la_factor_create(size_t m,size_t n,const double *input,int algorithm,la_factor **out) {
    if(!out||*out||algorithm<LA_LU||algorithm>LA_QR)return LA_ARGUMENT;
    if((algorithm!=LA_QR&&m!=n)||m<n)return LA_SHAPE;
    if((n&&m>SIZE_MAX/n/sizeof(double))||m>SIZE_MAX/sizeof(size_t)||n>SIZE_MAX/sizeof(double))return LA_MEMORY;
    if(m&&n&&!input)return LA_ARGUMENT;
    la_factor *f=(la_factor*)calloc(1,sizeof(la_factor));if(!f)return LA_MEMORY;
    f->rows=m;f->cols=n;f->algorithm=algorithm;f->scale=1;
    f->data=(double*)calloc(m*n?m*n:1,sizeof(double));f->tau=(double*)calloc(n?n:1,sizeof(double));
    f->permutation=(size_t*)calloc(m?m:1,sizeof(size_t));
    if(!f->data||!f->tau||!f->permutation){la_factor_destroy(f);return LA_MEMORY;}
    double largest=0;int status=LA_OK;double *a=f->data;
    for(size_t i=0;i<m*n;i++){if(!isfinite(input[i])){status=LA_RANGE;goto fail;}largest=fmax(largest,fabs(input[i]));}
    if(largest){int exponent;frexp(largest,&exponent);f->scale=ldexp(1,exponent-1);}
    for(size_t i=0;i<m*n;i++)if(la_scaled(input[i],f->scale,a+i)){status=LA_RANGE;goto fail;}
    for(size_t i=0;i<m;i++){
        double sum=0;for(size_t j=0;j<n;j++)sum+=fabs(a[i*n+j]);f->norm=fmax(f->norm,sum);f->permutation[i]=i;
    }
    if(algorithm==LA_LU){
        for(size_t k=0;k<n;k++){
            size_t pivot=k;for(size_t i=k+1;i<n;i++)if(fabs(a[i*n+k])>fabs(a[pivot*n+k]))pivot=i;
            if(a[pivot*n+k]==0){status=LA_SINGULAR;goto fail;}
            for(size_t j=0;j<n;j++){double t=a[k*n+j];a[k*n+j]=a[pivot*n+j];a[pivot*n+j]=t;}
            size_t p=f->permutation[k];f->permutation[k]=f->permutation[pivot];f->permutation[pivot]=p;
            for(size_t i=k+1;i<n;i++){
                a[i*n+k]/=a[k*n+k];if(!isfinite(a[i*n+k])){status=LA_RANGE;goto fail;}
                for(size_t j=k+1;j<n;j++){a[i*n+j]-=a[i*n+k]*a[k*n+j];if(!isfinite(a[i*n+j])){status=LA_RANGE;goto fail;}}
            }
        }
    }else if(algorithm==LA_CHOLESKY){
        for(size_t i=0;i<n;i++)for(size_t j=0;j<i;j++)if(input[i*n+j]!=input[j*n+i]){status=LA_NOT_SPD;goto fail;}
        for(size_t i=0;i<n;i++)for(size_t j=0;j<=i;j++){
            double value=a[i*n+j];for(size_t k=0;k<j;k++)value-=a[i*n+k]*a[j*n+k];
            if(!isfinite(value)){status=LA_RANGE;goto fail;}
            if(i==j){if(value<=0){status=LA_NOT_SPD;goto fail;}a[i*n+j]=sqrt(value);}
            else{a[i*n+j]=value/a[j*n+j];if(!isfinite(a[i*n+j])){status=LA_RANGE;goto fail;}}
        }
    }else{
        largest=0;for(size_t j=0;j<n;j++){double norm=0;for(size_t i=0;i<m;i++)norm=hypot(norm,a[i*n+j]);largest=fmax(largest,norm);}
        double threshold=DBL_EPSILON*(double)(m>n?m:n)*largest;
        for(size_t k=0;k<n;k++){
            size_t pivot=k;double norm=-1;
            for(size_t j=k;j<n;j++){double candidate=0;for(size_t i=k;i<m;i++)candidate=hypot(candidate,a[i*n+j]);if(candidate>norm){pivot=j;norm=candidate;}}
            if(norm<=threshold){status=LA_RANK;goto fail;}
            for(size_t i=0;i<m;i++){double t=a[i*n+k];a[i*n+k]=a[i*n+pivot];a[i*n+pivot]=t;}
            size_t p=f->permutation[k];f->permutation[k]=f->permutation[pivot];f->permutation[pivot]=p;
            double old=a[k*n+k],alpha=-copysign(norm,old),divisor=old-alpha;
            f->tau[k]=(alpha-old)/alpha;for(size_t i=k+1;i<m;i++)a[i*n+k]/=divisor;a[k*n+k]=alpha;
            for(size_t j=k+1;j<n;j++){
                double dot=a[k*n+j];for(size_t i=k+1;i<m;i++)dot+=a[i*n+k]*a[i*n+j];dot*=f->tau[k];
                a[k*n+j]-=dot;if(!isfinite(a[k*n+j])){status=LA_RANGE;goto fail;}
                for(size_t i=k+1;i<m;i++){a[i*n+j]-=a[i*n+k]*dot;if(!isfinite(a[i*n+j])){status=LA_RANGE;goto fail;}}
            }
        }
    }
    *out=f;return LA_OK;
fail:
    la_factor_destroy(f);return status;
}
/* Writes scratch output only. Public wrappers publish it after success. */
static int la_factor_solve(const la_factor *f,size_t rows,size_t p,const double *rhs,double *result,int rescale){
    if(!f)return LA_ARGUMENT;
    size_t m=f->rows,n=f->cols;if(rows!=m)return LA_SHAPE;
    if(p&&m>SIZE_MAX/p/sizeof(double))return LA_MEMORY;
    if((m&&p&&!rhs)||(n&&p&&!result))return LA_ARGUMENT;
    double *work=(double*)calloc(m*p?m*p:1,sizeof(double));if(!work)return LA_MEMORY;
    const double *a=f->data;int status=LA_OK;
    for(size_t i=0;i<m;i++)for(size_t j=0;j<p;j++){
        double value=rhs[(f->algorithm==LA_LU?f->permutation[i]:i)*p+j];
        if(!isfinite(value)){status=LA_RANGE;goto done;}
        if(rescale){status=la_scaled(value,f->scale,work+i*p+j);if(status)goto done;}else work[i*p+j]=value;
    }
    if(f->algorithm==LA_QR){
        for(size_t k=0;k<n;k++)for(size_t j=0;j<p;j++){
            double dot=work[k*p+j];for(size_t i=k+1;i<m;i++)dot+=a[i*n+k]*work[i*p+j];dot*=f->tau[k];
            work[k*p+j]-=dot;if(!isfinite(work[k*p+j])){status=LA_RANGE;goto done;}
            for(size_t i=k+1;i<m;i++){work[i*p+j]-=a[i*n+k]*dot;if(!isfinite(work[i*p+j])){status=LA_RANGE;goto done;}}
        }
    }else{
        for(size_t i=0;i<n;i++)for(size_t j=0;j<p;j++){
            double value=work[i*p+j];for(size_t k=0;k<i;k++)value-=a[i*n+k]*work[k*p+j];
            work[i*p+j]=f->algorithm==LA_CHOLESKY?value/a[i*n+i]:value;
            if(!isfinite(work[i*p+j])){status=LA_RANGE;goto done;}
        }
    }
    for(size_t i=n;i-->0;)for(size_t j=0;j<p;j++){
        double value=work[i*p+j];for(size_t k=i+1;k<n;k++)value-=(f->algorithm==LA_CHOLESKY?a[k*n+i]:a[i*n+k])*work[k*p+j];
        work[i*p+j]=value/a[i*n+i];if(!isfinite(work[i*p+j])){status=LA_RANGE;goto done;}
    }
    for(size_t i=0;i<n;i++)for(size_t j=0;j<p;j++)result[(f->algorithm==LA_QR?f->permutation[i]:i)*p+j]=work[i*p+j];
done:
    free(work);return status;
}
static int la_factor_rcond(const la_factor *f,double *out){
    if(!f||!out)return LA_ARGUMENT;
    size_t n=f->cols;if(f->rows!=n)return LA_SHAPE;if(!n){*out=1;return LA_OK;}
    double *identity=(double*)calloc(n*n,sizeof(double)),*inverse=(double*)calloc(n*n,sizeof(double));
    if(!identity||!inverse){free(identity);free(inverse);return LA_MEMORY;}
    for(size_t i=0;i<n;i++)identity[i*n+i]=1;
    int status=la_factor_solve(f,n,n,identity,inverse,0);double norm=0;
    if(status==LA_OK){for(size_t i=0;i<n;i++){double sum=0;for(size_t j=0;j<n;j++)sum+=fabs(inverse[i*n+j]);norm=fmax(norm,sum);}*out=fmin(1,(1/f->norm)/norm);}
    else if(status==LA_RANGE){*out=0;status=LA_OK;}
    free(identity);free(inverse);return status;
}
#endif
