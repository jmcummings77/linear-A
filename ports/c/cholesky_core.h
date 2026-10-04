#ifndef LINEAR_A_CHOLESKY_CORE_H
#define LINEAR_A_CHOLESKY_CORE_H
#include "sparse_core.h"
/* Owned plans and factors; initialize to zero, release with the matching free.
 * SIZE_MAX denotes original structure. Other steps are zero-based pivots. */
typedef struct {la_csr source,pattern;size_t *steps;size_t fill_count;} la_cholesky_symbolic;
typedef struct {la_csr lower;} la_cholesky;
typedef struct {size_t column,step;} la_cholesky_edge;
typedef struct {size_t count;la_cholesky_edge *edges;} la_cholesky_row;
static inline void la_cholesky_symbolic_free(la_cholesky_symbolic *s){if(s){la_csr_free(&s->source);la_csr_free(&s->pattern);free(s->steps);memset(s,0,sizeof(*s));}}
static inline void la_cholesky_free(la_cholesky *f){if(f)la_csr_free(&f->lower);}
static inline int la_cholesky_insert(la_cholesky_row *r,size_t j,size_t step){
 size_t k=0;while(k<r->count&&r->edges[k].column<j)k++;if(k<r->count&&r->edges[k].column==j)return 0;
 if(r->count==SIZE_MAX||!la_sparse_count(r->count+1,sizeof(la_cholesky_edge)))return 2;
 la_cholesky_edge *p=(la_cholesky_edge*)realloc(r->edges,(r->count+1)*sizeof(*p));if(!p)return 2;r->edges=p;
 memmove(p+k+1,p+k,(r->count-k)*sizeof(*p));p[k]=(la_cholesky_edge){j,step};r->count++;return 0;
}
static inline int la_cholesky_analyze_mode(const la_csr *a,la_cholesky_symbolic *out,int incomplete){
 if(!out||out->source.offsets||out->source.indices||out->source.values||out->source.rows||out->source.cols||out->source.nnz||out->pattern.offsets||out->pattern.indices||out->pattern.values||out->pattern.rows||out->pattern.cols||out->pattern.nnz||out->steps||out->fill_count||la_csr_validate(a)||a->rows!=a->cols)return 1;
 size_t n=a->rows;la_cholesky_symbolic s={0};int code=0;
 la_cholesky_row *g=(la_cholesky_row*)calloc(n?n:1,sizeof(*g));if(!g)return 2;
 for(size_t i=0;i<n&&!code;i++)for(size_t p=a->offsets[i];p<a->offsets[i+1]&&!code;p++){size_t j=a->indices[p];if(i!=j){code=la_cholesky_insert(&g[i],j,SIZE_MAX);if(!code)code=la_cholesky_insert(&g[j],i,SIZE_MAX);}}
 for(size_t k=0;!incomplete&&k<n&&!code;k++)for(size_t u=0;u<g[k].count&&!code;u++){size_t i=g[k].edges[u].column;if(i<=k)continue;for(size_t w=0;w<u&&!code;w++){size_t j=g[k].edges[w].column;if(j<=k)continue;code=la_cholesky_insert(&g[i],j,k);if(!code)code=la_cholesky_insert(&g[j],i,k);}}
 size_t count=n;
 for(size_t i=0;i<n&&!code;i++)for(size_t k=0;k<g[i].count;k++)if(g[i].edges[k].column<i){if(count==SIZE_MAX){code=2;break;}count++;}
 size_t *rp=NULL,*ci=NULL;double *v=NULL;
 if(!code){if(n==SIZE_MAX||!la_sparse_count(n+1,sizeof(size_t))||!la_sparse_count(count,sizeof(size_t))||!la_sparse_count(count,sizeof(double)))code=2;}
 if(!code){rp=(size_t*)calloc(n+1,sizeof(*rp));ci=(size_t*)calloc(count?count:1,sizeof(*ci));v=(double*)calloc(count?count:1,sizeof(*v));s.steps=(size_t*)calloc(count?count:1,sizeof(*s.steps));if(!rp||!ci||!v||!s.steps)code=2;}
 if(!code){size_t p=0;for(size_t i=0;i<n;i++){for(size_t k=0;k<g[i].count;k++){la_cholesky_edge e=g[i].edges[k];if(e.column<i){ci[p]=e.column;s.steps[p++]=e.step;if(e.step!=SIZE_MAX)s.fill_count++;}}ci[p]=i;s.steps[p++]=SIZE_MAX;rp[i+1]=p;}code=la_csr_create(n,n,count,rp,ci,v,&s.pattern);}
 /* Only the source pattern is retained; numerical values are not cached. */
 if(!code){double *zeros=(double*)calloc(a->nnz?a->nnz:1,sizeof(double));if(!zeros)code=2;else{code=la_csr_create(n,n,a->nnz,a->offsets,a->indices,zeros,&s.source);free(zeros);}}
 for(size_t i=0;i<n;i++)free(g[i].edges);free(g);free(rp);free(ci);free(v);
 if(code)la_cholesky_symbolic_free(&s);else *out=s;return code;
}
static inline int la_cholesky_analyze(const la_csr *a,la_cholesky_symbolic *out){return la_cholesky_analyze_mode(a,out,0);}
static inline int la_cholesky_factorize(const la_cholesky_symbolic *s,const la_csr *a,la_cholesky *out){
 if(!out||out->lower.offsets||out->lower.indices||out->lower.values||out->lower.rows||out->lower.cols||out->lower.nnz||!s||!s->steps||la_csr_validate(&s->source)||la_csr_validate(&s->pattern)||la_csr_validate(a)||a->rows!=a->cols||a->rows!=s->source.rows||a->nnz!=s->source.nnz||s->pattern.rows!=a->rows||s->pattern.cols!=a->rows)return 1;
 if(memcmp(a->offsets,s->source.offsets,(a->rows+1)*sizeof(size_t))||(a->nnz&&memcmp(a->indices,s->source.indices,a->nnz*sizeof(size_t))))return 1;
 for(size_t i=0;i<a->rows;i++){size_t end=s->pattern.offsets[i+1];if(end==s->pattern.offsets[i]||s->pattern.indices[end-1]!=i)return 1;for(size_t p=a->offsets[i];p<a->offsets[i+1];p++){size_t j=a->indices[p],q=la_csr_find(a,j,i);double t=q<a->offsets[j+1]&&a->indices[q]==i?a->values[q]:0;if(a->values[p]!=t)return 4;}}
 la_cholesky f={0};int code=la_csr_create(a->rows,a->cols,s->pattern.nnz,s->pattern.offsets,s->pattern.indices,s->pattern.values,&f.lower);if(code)return code;la_csr *l=&f.lower;
 for(size_t i=0;i<a->rows&&!code;i++)for(size_t p=l->offsets[i];p<l->offsets[i+1];p++){
  size_t j=l->indices[p],q=la_csr_find(a,i,j),u=l->offsets[i],w=l->offsets[j];double v=q<a->offsets[i+1]&&a->indices[q]==j?a->values[q]:0;
  while(u<p&&w<l->offsets[j+1]-1){if(l->indices[u]==l->indices[w]){v-=l->values[u]*l->values[w];u++;w++;}else if(l->indices[u]<l->indices[w])u++;else w++;}
  if(!isfinite(v)){code=3;break;}if(i==j){if(v<=0){code=5;break;}l->values[p]=sqrt(v);}else{l->values[p]=v/l->values[l->offsets[j+1]-1];if(!isfinite(l->values[p])){code=3;break;}}
 }
 if(code)la_cholesky_free(&f);else *out=f;return code;
}
static inline int la_cholesky_solve(const la_cholesky *f,const double *b,size_t count,double *x){
 if(!f||la_csr_validate(&f->lower)||f->lower.rows!=f->lower.cols||count!=f->lower.rows||(!b&&count)||(!x&&count))return 1;const la_csr *l=&f->lower;
 for(size_t i=0;i<count;i++){size_t end=l->offsets[i+1];if(!isfinite(b[i])||end==l->offsets[i]||l->indices[end-1]!=i||l->values[end-1]<=0)return 1;}
 if(count)memmove(x,b,count*sizeof(double));
 for(size_t i=0;i<count;i++){size_t d=l->offsets[i+1]-1;for(size_t p=l->offsets[i];p<d;p++)x[i]-=l->values[p]*x[l->indices[p]];x[i]/=l->values[d];if(!isfinite(x[i]))return 3;}
 for(size_t ii=count;ii>0;ii--){size_t i=ii-1,d=l->offsets[i+1]-1;x[i]/=l->values[d];if(!isfinite(x[i]))return 3;for(size_t p=l->offsets[i];p<d;p++){x[l->indices[p]]-=l->values[p]*x[i];if(!isfinite(x[l->indices[p]]))return 3;}}
 return 0;
}
typedef la_cholesky la_ic0;
static inline int la_ic0_create(const la_csr *a,la_ic0 *out){la_cholesky_symbolic s={0};int code=la_cholesky_analyze_mode(a,&s,1);if(!code)code=la_cholesky_factorize(&s,a,out);la_cholesky_symbolic_free(&s);return code;}
static inline int la_ic0_apply(const la_ic0 *f,const double *b,size_t n,double *x){return la_cholesky_solve(f,b,n,x);}
static inline int la_ic0_callback(const void *f,const double *b,size_t n,double *x){return la_ic0_apply((const la_ic0*)f,b,n,x);}
static inline int la_csr_cg_ic0(const la_csr *a,const double *b,size_t count,double rtol,double atol,size_t limit,int jacobi,int capture,const la_ic0 *f,la_cg_result *out){
 if(f){if(!a||jacobi||la_csr_validate(&f->lower)||f->lower.rows!=a->rows||f->lower.cols!=a->rows)return 1;
 for(size_t i=0;i<f->lower.rows;i++){size_t end=f->lower.offsets[i+1];if(end==f->lower.offsets[i]||f->lower.indices[end-1]!=i||f->lower.values[end-1]<=0)return 1;}}
 return la_csr_cg_with_apply(a,b,count,rtol,atol,limit,jacobi,capture,f?la_ic0_callback:NULL,f,out);
}
#endif
