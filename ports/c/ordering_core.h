#ifndef LINEAR_A_ORDERING_CORE_H
#define LINEAR_A_ORDERING_CORE_H
#include "sparse_core.h"
typedef struct {size_t first,second;} la_order_pair;
static int la_order_compare(const void *aa,const void *bb){const la_order_pair *a=(const la_order_pair*)aa,*b=(const la_order_pair*)bb;return a->first<b->first?-1:a->first>b->first?1:a->second<b->second?-1:a->second>b->second;}
/* Minimum-degree component seed, degree/index neighbor ties; reverse entire BFS order. */
static inline int la_csr_rcm(const la_csr *a,size_t *out,size_t count){
 if(la_csr_validate(a)||a->rows!=a->cols||count!=a->rows||(!out&&count))return 1;
 size_t n=a->rows;
 if(!la_sparse_count(a->nnz,2*sizeof(la_order_pair))||!la_sparse_count(n,sizeof(la_order_pair))||!la_sparse_count(n+1,sizeof(size_t)))return 2;
 la_order_pair *edges=(la_order_pair*)calloc(a->nnz?2*a->nnz:1,sizeof(*edges)),*starts=(la_order_pair*)calloc(n?n:1,sizeof(*starts)),*next=(la_order_pair*)calloc(n?n:1,sizeof(*next));
 size_t *rp=(size_t*)calloc(n+1,sizeof(size_t)),*queue=(size_t*)calloc(n?n:1,sizeof(size_t));unsigned char *seen=(unsigned char*)calloc(n?n:1,1);
 if(!edges||!starts||!next||!rp||!queue||!seen){free(edges);free(starts);free(next);free(rp);free(queue);free(seen);return 2;}
 size_t m=0;
 for(size_t i=0;i<n;i++)for(size_t k=a->offsets[i];k<a->offsets[i+1];k++){size_t j=a->indices[k];if(i!=j){edges[m++]=(la_order_pair){i,j};edges[m++]=(la_order_pair){j,i};}}
 qsort(edges,m,sizeof(*edges),la_order_compare);size_t unique=0;
 for(size_t k=0;k<m;k++)if(!unique||la_order_compare(edges+k,edges+unique-1)){edges[unique++]=edges[k];}
 for(size_t k=0;k<unique;k++)rp[edges[k].first+1]++;
 for(size_t i=0;i<n;i++){rp[i+1]+=rp[i];starts[i]=(la_order_pair){rp[i+1]-rp[i],i};}
 qsort(starts,n,sizeof(*starts),la_order_compare);size_t used=0;
 for(size_t s=0;s<n;s++){size_t start=starts[s].second;if(seen[start])continue;size_t head=used;queue[used++]=start;seen[start]=1;
  while(head<used){size_t i=queue[head++],length=0;
   for(size_t k=rp[i];k<rp[i+1];k++){size_t j=edges[k].second;if(!seen[j])next[length++]=(la_order_pair){rp[j+1]-rp[j],j};}
   qsort(next,length,sizeof(*next),la_order_compare);for(size_t k=0;k<length;k++){size_t j=next[k].second;seen[j]=1;queue[used++]=j;}
  }
 }
 for(size_t i=0;i<n;i++)out[i]=queue[n-1-i];
 free(edges);free(starts);free(next);free(rp);free(queue);free(seen);return 0;
}
static inline int la_permutation_inverse(const size_t *p,size_t n,size_t *inv){
 if((!p||!inv)&&n)return 1;
 for(size_t i=0;i<n;i++)inv[i]=SIZE_MAX;
 for(size_t i=0;i<n;i++){if(p[i]>=n||inv[p[i]]!=SIZE_MAX)return 1;inv[p[i]]=i;}return 0;
}
static inline int la_permute_vector(const size_t *p,size_t n,const double *x,int inverse,double *out){
 if(((!p||!x||!out)&&n)||(inverse!=0&&inverse!=1))return 1;
 if(!la_sparse_count(n,sizeof(size_t))||!la_sparse_count(n,sizeof(double)))return 2;
 size_t *inv=(size_t*)calloc(n?n:1,sizeof(size_t));double *tmp=(double*)calloc(n?n:1,sizeof(double));if(!inv||!tmp){free(inv);free(tmp);return 2;}
 int code=la_permutation_inverse(p,n,inv);
 for(size_t i=0;!code&&i<n;i++)if(!isfinite(x[i]))code=1;
 if(!code){for(size_t i=0;i<n;i++)tmp[i]=x[inverse?inv[i]:p[i]];if(n)memcpy(out,tmp,n*sizeof(double));}
 free(inv);free(tmp);return code;
}
static inline int la_csr_permute(const la_csr *a,const size_t *p,size_t count,la_csr *out){
 if(la_csr_validate(a)||a->rows!=a->cols||count!=a->rows||!out||out->offsets||out->indices||out->values||out->rows||out->cols||out->nnz)return 1;
 size_t n=a->rows;if(!la_sparse_count(n,sizeof(size_t))||!la_sparse_count(a->nnz,sizeof(la_order_pair)))return 2;
 size_t *inv=(size_t*)calloc(n?n:1,sizeof(size_t));la_order_pair *entries=(la_order_pair*)calloc(a->nnz?a->nnz:1,sizeof(*entries));if(!inv||!entries){free(inv);free(entries);return 2;}
 int code=la_permutation_inverse(p,n,inv);la_csr result={0};if(!code)code=la_csr_create(n,n,a->nnz,a->offsets,a->indices,a->values,&result);
 if(!code){size_t used=0;result.offsets[0]=0;
  for(size_t row=0;row<n;row++){size_t i=p[row],length=0;for(size_t k=a->offsets[i];k<a->offsets[i+1];k++)entries[length++]=(la_order_pair){inv[a->indices[k]],k};qsort(entries,length,sizeof(*entries),la_order_compare);
   for(size_t k=0;k<length;k++){result.indices[used]=entries[k].first;result.values[used++]=a->values[entries[k].second];}result.offsets[row+1]=used;
  }*out=result;
 }
 free(inv);free(entries);return code;
}
#endif
