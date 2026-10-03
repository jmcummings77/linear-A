#ifndef LINEAR_A_SVD_CORE_H
#define LINEAR_A_SVD_CORE_H
/* Economy one-sided Jacobi SVD. Shared by C, C++, ARM64 and WebAssembly.
 * Input/output storage is row-major; outputs must be distinct allocated arrays.
 * Status: 0 success, 1 argument, 2 allocation, 3 range, 4 nonconvergence. */
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
static double la_svd_norm(const double *b,size_t m,size_t n,size_t j){
    double norm=0;for(size_t i=0;i<m;i++)norm=hypot(norm,b[i*n+j]);return norm;
}
static int la_svd(size_t m,size_t n,const double *input,double tolerance,size_t sweeps,double *u,double *values,double *vt){
    if(!isfinite(tolerance)||tolerance<=0||tolerance>=1||!sweeps||sweeps>10000)return 1;
    if(!m||!n)return 0;
    if(m<n){
        if(m&&n>SIZE_MAX/m/sizeof(double))return 2;
        double *at=(double*)calloc(m*n?m*n:1,sizeof(double));
        double *tu=(double*)calloc(m*n?m*n:1,sizeof(double));
        double *tv=(double*)calloc(m*m?m*m:1,sizeof(double));
        if(!at||!tu||!tv){free(at);free(tu);free(tv);return 2;}
        for(size_t i=0;i<m;i++)for(size_t j=0;j<n;j++)at[j*m+i]=input[i*n+j];
        int status=la_svd(n,m,at,tolerance,sweeps,tu,values,tv);
        if(!status){
            for(size_t i=0;i<m;i++)for(size_t j=0;j<m;j++)u[i*m+j]=tv[j*m+i];
            for(size_t i=0;i<m;i++)for(size_t j=0;j<n;j++)vt[i*n+j]=tu[j*m+i];
        }
        free(at);free(tu);free(tv);return status;
    }
    if(n&&m>SIZE_MAX/n/sizeof(double))return 2;
    double *b=(double*)calloc(m*n?m*n:1,sizeof(double)),*v=(double*)calloc(n*n?n*n:1,sizeof(double));
    double *norms=(double*)calloc(n?n:1,sizeof(double)),*candidate=(double*)calloc(m?m:1,sizeof(double));
    size_t *order=(size_t*)calloc(n?n:1,sizeof(size_t));
    if(!b||!v||!norms||!candidate||!order){free(b);free(v);free(norms);free(candidate);free(order);return 2;}
    int status=0,converged=0;double scale=0;
    for(size_t i=0;i<m*n;i++){if(!isfinite(input[i])){status=3;goto done;}scale=fmax(scale,fabs(input[i]));}
    if(!scale)scale=1;
    for(size_t i=0;i<m*n;i++){b[i]=input[i]/scale;if(input[i]&&!b[i]){status=3;goto done;}}
    for(size_t j=0;j<n;j++)v[j*n+j]=1;
    for(size_t sweep=0;sweep<=sweeps;sweep++){
        int changed=0;
        for(size_t p=0;p<n;p++)for(size_t q=p+1;q<n;q++){
            double np=la_svd_norm(b,m,n,p),nq=la_svd_norm(b,m,n,q),corr=0;
            if(!np||!nq)continue;
            for(size_t i=0;i<m;i++)corr+=(b[i*n+p]/np)*(b[i*n+q]/nq);
            if(fabs(corr)<=tolerance)continue;
            changed=1;if(sweep==sweeps)continue;
            double pair=fmax(np,nq),ap=np/pair,aq=nq/pair,delta=aq*aq-ap*ap,g=2*ap*aq*corr;
            double t=delta==0?copysign(1,g):g/(delta+copysign(hypot(delta,g),delta));
            if(fabs(t)<DBL_MIN){size_t small=np<nq?p:q;for(size_t i=0;i<m;i++)b[i*n+small]=0;continue;}
            double c=1/hypot(1,t),s=c*t;
            for(size_t i=0;i<m;i++){double x=b[i*n+p],y=b[i*n+q];b[i*n+p]=c*x-s*y;b[i*n+q]=s*x+c*y;}
            for(size_t i=0;i<n;i++){double x=v[i*n+p],y=v[i*n+q];v[i*n+p]=c*x-s*y;v[i*n+q]=s*x+c*y;}
        }
        if(!changed){converged=1;break;}
    }
    if(!converged){status=4;goto done;}
    for(size_t j=0;j<n;j++){norms[j]=la_svd_norm(b,m,n,j);order[j]=j;}
    for(size_t j=1;j<n;j++){size_t x=order[j],k=j;while(k&&norms[order[k-1]]<norms[x]){order[k]=order[k-1];k--;}order[k]=x;}
    for(size_t j=0;j<n;j++){
        size_t k=order[j];values[j]=norms[k]*scale;
        if(!isfinite(values[j])||(norms[k]&&!values[j])){status=3;goto done;}
        for(size_t i=0;i<n;i++)vt[j*n+i]=v[i*n+k];
        if(norms[k]){for(size_t i=0;i<m;i++)u[i*n+j]=b[i*n+k]/norms[k];}
        else{
            int found=0;
            for(size_t axis=0;axis<m;axis++){
                for(size_t i=0;i<m;i++)candidate[i]=(double)(i==axis);
                for(int pass=0;pass<2;pass++)for(size_t col=0;col<j;col++){
                    double dot=0;for(size_t i=0;i<m;i++)dot+=candidate[i]*u[i*n+col];
                    for(size_t i=0;i<m;i++)candidate[i]-=dot*u[i*n+col];
                }
                double length=0;for(size_t i=0;i<m;i++)length=hypot(length,candidate[i]);
                if(length>0.5/sqrt((double)m)){for(size_t i=0;i<m;i++)u[i*n+j]=candidate[i]/length;found=1;break;}
            }
            if(!found){status=4;goto done;}
        }
    }
 done:free(b);free(v);free(norms);free(candidate);free(order);return status;
}
#endif
