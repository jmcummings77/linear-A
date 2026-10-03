#ifndef LINEAR_A_GMRES_CORE_H
#define LINEAR_A_GMRES_CORE_H
#include "sparse_core.h"
#include <float.h>
typedef struct { double *x,*residuals,*estimated_residuals,*iterates; size_t *restarts; size_t size,iterations,restart_count; int reason; } la_gmres_result;
static inline void la_gmres_free(la_gmres_result *r){if(r){free(r->x);free(r->residuals);free(r->estimated_residuals);free(r->iterates);free(r->restarts);memset(r,0,sizeof(*r));}}
static inline int la_gmres_finite(size_t n,const double *v){for(size_t i=0;i<n;i++)if(!isfinite(v[i]))return 0;return 1;}
/* Right Jacobi preconditioning keeps projected and true residuals in the same norm.
 * Reason codes: 0 converged, 1 limit, 2 breakdown, 3 nonfinite, 4 stagnation. */
static inline int la_csr_gmres(const la_csr *a,const double *b,size_t count,size_t restart,double rtol,double atol,size_t limit,int jacobi,int capture,la_gmres_result *out){
    if(!out||out->x||out->residuals||out->estimated_residuals||out->iterates||out->restarts||out->size||out->iterations||out->restart_count||out->reason)return 1;
    int code=la_csr_validate(a);if(code)return code;
    size_t n=a->rows,m=restart;if(m>n)m=n;if(m>limit)m=limit;
    if(n!=a->cols||count!=n||(!b&&n)||restart<1||restart>1024||limit>100000||!isfinite(rtol)||rtol<0||rtol>=1||!isfinite(atol)||atol<0||(jacobi!=0&&jacobi!=1)||(capture!=0&&capture!=1)||!la_gmres_finite(n,b))return 1;
    if(!la_sparse_count(n,8*sizeof(double))||n>SIZE_MAX/(m+1)/sizeof(double)||(capture&&n>SIZE_MAX/(limit+1)/sizeof(double)))return 2;
    double *work=(double*)calloc(n?n*8:1,sizeof(double)),*basis=(double*)calloc(n?(m+1)*n:1,sizeof(double)),*h=(double*)calloc((m+1)*(m?m:1),sizeof(double));
    double *cs=(double*)calloc(m+1,sizeof(double)),*sn=(double*)calloc(m+1,sizeof(double)),*g=(double*)calloc(m+1,sizeof(double)),*y=(double*)calloc(m+1,sizeof(double));
    if(!work||!basis||!h||!cs||!sn||!g||!y){free(work);free(basis);free(h);free(cs);free(sn);free(g);free(y);return 2;}
    double *r=work,*base=r+n,*w=base+n,*z=w+n,*candidate=z+n,*ax=candidate+n,*res=ax+n,*diag=res+n;
    for(size_t i=0;i<n;i++){
        diag[i]=1;
        if(jacobi){size_t p=la_csr_find(a,i,i);if(p==a->offsets[i+1]||a->indices[p]!=i||a->values[p]==0){code=1;goto cleanup;}diag[i]=a->values[p];}
    }
    {
        la_gmres_result result={0};result.size=n;result.reason=1;
        result.x=(double*)calloc(n?n:1,sizeof(double));result.residuals=(double*)calloc(limit+1,sizeof(double));result.estimated_residuals=(double*)calloc(limit+1,sizeof(double));result.restarts=(size_t*)calloc(limit+1,sizeof(size_t));
        if(capture)result.iterates=(double*)calloc(n?(limit+1)*n:1,sizeof(double));
        if(!result.x||!result.residuals||!result.estimated_residuals||!result.restarts||(capture&&!result.iterates)){la_gmres_free(&result);code=2;goto cleanup;}
        if(n)memcpy(r,b,n*sizeof(double));result.residuals[0]=result.estimated_residuals[0]=la_sparse_norm(n,r);
        double threshold=fmax(atol,rtol*result.residuals[0]);
        if(!isfinite(result.residuals[0])){result.reason=3;goto done;}
        if(result.residuals[0]<=threshold){result.reason=0;goto done;}
        while(result.iterations<limit){
            if(result.iterations)result.restarts[result.restart_count++]=result.iterations;
            memcpy(base,result.x,n*sizeof(double));double beta=la_sparse_norm(n,r);
            memset(h,0,(m+1)*m*sizeof(double));memset(g,0,(m+1)*sizeof(double));g[0]=beta;
            for(size_t i=0;i<n;i++)basis[i]=r[i]/beta;
            size_t steps=m;if(steps>limit-result.iterations)steps=limit-result.iterations;
            for(size_t j=0;j<steps;j++){
                for(size_t i=0;i<n;i++)z[i]=basis[j*n+i]/diag[i];
                if(la_csr_mv(a,z,w)){result.reason=3;goto done;}
                double original=la_sparse_norm(n,w);
                for(int pass=0;pass<2;pass++)for(size_t k=0;k<=j;k++){
                    double dot=la_sparse_dot(n,basis+k*n,w);h[k*m+j]+=dot;
                    for(size_t i=0;i<n;i++)w[i]-=dot*basis[k*n+i];
                }
                double tail=la_sparse_norm(n,w);
                if(!isfinite(original)||!isfinite(tail)){result.reason=3;goto done;}
                for(size_t k=0;k<=j;k++)if(!isfinite(h[k*m+j])){result.reason=3;goto done;}
                int happy=tail<=8*DBL_EPSILON*original;h[(j+1)*m+j]=happy?0:tail;
                if(!happy)for(size_t i=0;i<n;i++)basis[(j+1)*n+i]=w[i]/tail;
                for(size_t k=0;k<j;k++){
                    double top=cs[k]*h[k*m+j]+sn[k]*h[(k+1)*m+j];
                    h[(k+1)*m+j]=-sn[k]*h[k*m+j]+cs[k]*h[(k+1)*m+j];h[k*m+j]=top;
                }
                double pivot=hypot(h[j*m+j],h[(j+1)*m+j]);
                if(!isfinite(pivot)){result.reason=3;goto done;}if(pivot==0){result.reason=2;goto done;}
                cs[j]=h[j*m+j]/pivot;sn[j]=h[(j+1)*m+j]/pivot;h[j*m+j]=pivot;h[(j+1)*m+j]=0;
                g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];memcpy(y,g,(j+1)*sizeof(double));
                for(size_t kk=j+1;kk>0;kk--){size_t k=kk-1;if(h[k*m+k]==0){result.reason=2;goto done;}double sum=0;for(size_t q=k+1;q<=j;q++)sum+=h[k*m+q]*y[q];y[k]=(y[k]-sum)/h[k*m+k];}
                for(size_t i=0;i<n;i++){double sum=0;for(size_t k=0;k<=j;k++)sum+=basis[k*n+i]*y[k];candidate[i]=base[i]+sum/diag[i];}
                if(!la_gmres_finite(j+1,y)||!isfinite(g[j+1])||la_csr_mv(a,candidate,ax)){result.reason=3;goto done;}
                for(size_t i=0;i<n;i++)res[i]=b[i]-ax[i];double length=la_sparse_norm(n,res);
                if(!isfinite(length)){result.reason=3;goto done;}
                memcpy(result.x,candidate,n*sizeof(double));memcpy(r,res,n*sizeof(double));result.iterations++;
                result.residuals[result.iterations]=length;result.estimated_residuals[result.iterations]=fabs(g[j+1]);
                if(capture)memcpy(result.iterates+result.iterations*n,result.x,n*sizeof(double));
                if(length<=threshold){result.reason=0;goto done;}if(happy){result.reason=2;goto done;}
            }
            {int same=1;for(size_t i=0;i<n;i++)if(result.x[i]!=base[i])same=0;if(same){result.reason=4;goto done;}}
        }
    done:*out=result;
    }
cleanup:free(work);free(basis);free(h);free(cs);free(sn);free(g);free(y);return code;
}
#endif
