#ifndef LINEAR_A_MULTIGRID_CORE_H
#define LINEAR_A_MULTIGRID_CORE_H
#include "sparse_core.h"
/* Unit five-point Dirichlet Laplacian; geometric nested grids. */
typedef struct {size_t width;} la_multigrid;
typedef int (*la_mg_record)(void*,size_t,size_t,int,const double*,const double*,const double*);
static inline int la_mg_valid(const la_multigrid *m){return m&&m->width>=1&&m->width<=255&&!(m->width&(m->width+1));}
static inline int la_mg_create(size_t width,la_multigrid *out){la_multigrid m={width};if(!out||out->width||!la_mg_valid(&m))return 1;*out=m;return 0;}
static inline void la_mg_mv(size_t w,const double *x,double *out){for(size_t i=0;i<w*w;i++){size_t y=i/w,j=i%w;out[i]=4*x[i]-(j?x[i-1]:0)-(j+1<w?x[i+1]:0)-(y?x[i-w]:0)-(y+1<w?x[i+w]:0);}}
static inline int la_mg_emit(size_t level,size_t w,int phase,const double *x,const double *b,double *r,la_mg_record record,void *context){size_t n=w*w;for(size_t i=0;i<n;i++)if(!isfinite(x[i]))return 3;if(!record)return 0;la_mg_mv(w,x,r);for(size_t i=0;i<n;i++){r[i]=b[i]-r[i];if(!isfinite(r[i]))return 3;}return record(context,level,w,phase,x,b,r);}
static inline void la_mg_smooth(size_t w,const double *b,double *x,double *tmp){for(int k=0;k<2;k++){la_mg_mv(w,x,tmp);for(size_t i=0;i<w*w;i++)x[i]+=(b[i]-tmp[i])/6;}}
static inline int la_mg_cycle(size_t w,const double *b,double *x,size_t level,la_mg_record record,void *context){
 size_t n=w*w,c=w/2,nc=c*c;double *work=(double*)calloc(n+2*nc,sizeof(double));if(!work)return 2;double *r=work,*bc=r+n,*ec=bc+nc;memset(x,0,n*sizeof(double));int code=la_mg_emit(level,w,0,x,b,r,record,context);
 if(!code&&w==1){x[0]=b[0]/4;code=la_mg_emit(level,w,2,x,b,r,record,context);free(work);return code;}
 if(!code){la_mg_smooth(w,b,x,r);code=la_mg_emit(level,w,1,x,b,r,record,context);}
 if(!code){la_mg_mv(w,x,r);for(size_t i=0;i<n;i++){r[i]=b[i]-r[i];if(!isfinite(r[i]))code=3;}}
 if(!code){for(size_t y=0;y<c;y++)for(size_t j=0;j<c;j++)for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++)bc[y*c+j]+=(dy==0?1:.5)*(dx==0?1:.5)*r[(2*y+1+dy)*w+2*j+1+dx];for(size_t i=0;i<nc;i++)if(!isfinite(bc[i]))code=3;}
 if(!code)code=la_mg_cycle(c,bc,ec,level+1,record,context);
 if(!code){for(size_t y=0;y<c;y++)for(size_t j=0;j<c;j++)for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++)x[(2*y+1+dy)*w+2*j+1+dx]+=(dy==0?1:.5)*(dx==0?1:.5)*ec[y*c+j];code=la_mg_emit(level,w,3,x,b,r,record,context);}
 if(!code){la_mg_smooth(w,b,x,r);code=la_mg_emit(level,w,4,x,b,r,record,context);}free(work);return code;
}
static inline int la_mg_apply_trace(const la_multigrid *m,const double *b,size_t count,double *x,la_mg_record record,void *context){if(!la_mg_valid(m)||count!=m->width*m->width||!b||!x)return 1;for(size_t i=0;i<count;i++)if(!isfinite(b[i]))return 1;double *result=(double*)calloc(count,sizeof(double));if(!result)return 2;int code=la_mg_cycle(m->width,b,result,0,record,context);if(!code)memcpy(x,result,count*sizeof(double));free(result);return code;}
static inline int la_mg_apply(const la_multigrid *m,const double *b,size_t count,double *x){return la_mg_apply_trace(m,b,count,x,NULL,NULL);}
static inline int la_mg_callback(const void *m,const double *b,size_t n,double *x){return la_mg_apply((const la_multigrid*)m,b,n,x);}
static inline int la_csr_cg_mg(const la_csr *a,const double *b,size_t count,double rtol,double atol,size_t limit,int jacobi,int capture,const la_multigrid *m,la_cg_result *out){if(!a||!la_mg_valid(m)||a->rows!=m->width*m->width||jacobi)return 1;return la_csr_cg_with_apply(a,b,count,rtol,atol,limit,jacobi,capture,la_mg_callback,m,out);}
static inline int la_mg_matrix(const la_multigrid *m,la_csr *out){if(!la_mg_valid(m))return 1;size_t w=m->width,n=w*w,nnz=5*n-4*w;size_t *rp=(size_t*)calloc(n+1,sizeof(size_t)),*ci=(size_t*)calloc(nnz,sizeof(size_t));double *v=(double*)calloc(nnz,sizeof(double));if(!rp||!ci||!v){free(rp);free(ci);free(v);return 2;}size_t p=0;for(size_t i=0;i<n;i++){if(i>=w){ci[p]=i-w;v[p++]=-1;}if(i%w){ci[p]=i-1;v[p++]=-1;}ci[p]=i;v[p++]=4;if(i%w+1<w){ci[p]=i+1;v[p++]=-1;}if(i+w<n){ci[p]=i+w;v[p++]=-1;}rp[i+1]=p;}int code=la_csr_create(n,n,nnz,rp,ci,v,out);free(rp);free(ci);free(v);return code;}
#endif
