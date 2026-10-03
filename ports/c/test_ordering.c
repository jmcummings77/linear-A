#include "matrix.h"
#include <assert.h>
#include <math.h>
int main(void){
 size_t rp[]={0,2,3,5},ci[]={0,2,1,0,2},p[]={2,0,1};double v[]={4,0,5,2,6},x[]={1,2,3},y[3],back[3];sparse_matrix a={0},q={0};
 assert(m_csr_create(3,3,5,rp,ci,v,&a)==M_OK);assert(m_csr_permute(&a,p,3,&q)==M_OK);assert(q.nnz==5&&q.values[1]==2&&q.values[2]==0);
 assert(m_permute_vector(p,3,x,false,y)==M_OK);assert(y[0]==3&&y[1]==1&&y[2]==2);assert(m_permute_vector(p,3,y,true,back)==M_OK);for(int i=0;i<3;i++)assert(x[i]==back[i]);
 assert(m_permute_vector(p,3,x,false,x)==M_OK);assert(x[0]==3&&x[1]==1&&x[2]==2);
 size_t bad[]={0,0,2};sparse_matrix invalid={0};assert(m_csr_permute(&a,bad,3,&invalid)==M_ARGUMENT&&!invalid.offsets);double untouched[]={9,9,9};assert(m_permute_vector(bad,3,x,false,untouched)==M_ARGUMENT&&untouched[0]==9);
 assert(m_csr_permute(&a,p,2,&invalid)==M_ARGUMENT);assert(m_csr_permute(&a,p,3,&q)==M_ARGUMENT);
 size_t order[3];assert(m_csr_rcm(&a,order,3)==M_OK);assert(order[0]==2&&order[1]==0&&order[2]==1);assert(m_csr_rcm(&a,order,2)==M_ARGUMENT);
 m_csr_free(&a);assert(q.values[0]==6);m_csr_free(&q);m_csr_free(&q);
 size_t empty[]={0};assert(m_csr_create(0,0,0,empty,NULL,NULL,&a)==M_OK);assert(m_csr_rcm(&a,NULL,0)==M_OK);assert(m_csr_permute(&a,NULL,0,&q)==M_OK);m_csr_free(&a);m_csr_free(&q);return 0;
}
