#ifndef LINEAR_A_ILU_CORE_H
#define LINEAR_A_ILU_CORE_H
#include "sparse_core.h"
/* Owning snapshot. Zero fill, no permutations or shifts. Code 6: zero pivot. */
typedef struct { la_csr factors; size_t *diagonal; } la_ilu0;
static inline void la_ilu0_free(la_ilu0 *f){if(f){la_csr_free(&f->factors);free(f->diagonal);memset(f,0,sizeof(*f));}}
static inline int la_ilu0_validate(const la_ilu0 *f){
 if(!f||la_csr_validate(&f->factors)||f->factors.rows!=f->factors.cols||!f->diagonal)return 1;
 const la_csr *a=&f->factors;
 for(size_t i=0;i<a->rows;i++){size_t p=f->diagonal[i];if(p<a->offsets[i]||p>=a->offsets[i+1]||a->indices[p]!=i||a->values[p]==0)return 1;}return 0;
}
static inline int la_ilu0_create(const la_csr *source,la_ilu0 *out){
 if(!out||out->diagonal||out->factors.offsets||out->factors.indices||out->factors.values||out->factors.rows||out->factors.cols||out->factors.nnz||la_csr_validate(source)||source->rows!=source->cols)return 1;
 la_ilu0 f={0};int code=la_csr_create(source->rows,source->cols,source->nnz,source->offsets,source->indices,source->values,&f.factors);if(code)return code;
 la_csr *a=&f.factors;size_t n=a->rows;f.diagonal=(size_t*)calloc(n?n:1,sizeof(size_t));if(!f.diagonal){la_ilu0_free(&f);return 2;}
 for(size_t i=0;i<n;i++){size_t p=la_csr_find(a,i,i);if(p==a->offsets[i+1]||a->indices[p]!=i){la_ilu0_free(&f);return 1;}f.diagonal[i]=p;}
 for(size_t i=0;i<n;i++){
  for(size_t p=a->offsets[i];p<f.diagonal[i];p++){size_t j=a->indices[p];a->values[p]/=a->values[f.diagonal[j]];
   if(!isfinite(a->values[p])){la_ilu0_free(&f);return 3;}
   for(size_t q=f.diagonal[j]+1;q<a->offsets[j+1];q++){size_t k=la_csr_find(a,i,a->indices[q]);if(k<a->offsets[i+1]&&a->indices[k]==a->indices[q]){a->values[k]-=a->values[p]*a->values[q];if(!isfinite(a->values[k])){la_ilu0_free(&f);return 3;}}}
  }
  if(a->values[f.diagonal[i]]==0){la_ilu0_free(&f);return 6;}
 }
 *out=f;return 0;
}
static inline int la_ilu0_apply(const la_ilu0 *f,const double *b,size_t count,double *x){
 if(la_ilu0_validate(f)||count!=f->factors.rows||(!b&&count)||(!x&&count))return 1;
 const la_csr *a=&f->factors;size_t n=a->rows;
 for(size_t i=0;i<n;i++)if(!isfinite(b[i]))return 1;
 if(n)memmove(x,b,n*sizeof(double));
 for(size_t i=0;i<n;i++){for(size_t p=a->offsets[i];p<f->diagonal[i];p++)x[i]-=a->values[p]*x[a->indices[p]];if(!isfinite(x[i]))return 3;}
 for(size_t ii=n;ii>0;ii--){size_t i=ii-1;for(size_t p=f->diagonal[i]+1;p<a->offsets[i+1];p++)x[i]-=a->values[p]*x[a->indices[p]];x[i]/=a->values[f->diagonal[i]];if(!isfinite(x[i]))return 3;}
 return 0;
}
#endif
