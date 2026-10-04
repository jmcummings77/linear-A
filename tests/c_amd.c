#include "../ports/c/matrix.h"
#include <assert.h>
int main(void){
 size_t rp[]={0,3,5,7},ci[]={0,1,2,0,1,0,2},p[3];double v[]={4,-1,-1,-1,2,-1,2};sparse_matrix a={0};
 assert(m_csr_create(3,3,7,rp,ci,v,&a)==M_OK);
 assert(m_csr_amd(&a,p,3)==M_OK&&p[0]==1&&p[1]==0&&p[2]==2);
 assert(m_csr_amd(&a,p,2)!=M_OK);assert(m_csr_amd(&a,NULL,3)!=M_OK);
 a.cols=4;assert(m_csr_amd(&a,p,3)!=M_OK);a.cols=3;
 m_csr_free(&a);size_t empty[]={0};assert(m_csr_create(0,0,0,empty,NULL,NULL,&a)==M_OK);assert(m_csr_amd(&a,NULL,0)==M_OK);m_csr_free(&a);
 /* Dense graph exercises set growth, clique absorption and deterministic ties. */
 size_t offsets[21],indices[400],order[20];double values[400];
 for(size_t i=0;i<20;i++){offsets[i]=20*i;for(size_t j=0;j<20;j++){indices[20*i+j]=j;values[20*i+j]=i==j?21:-1;}}offsets[20]=400;
 assert(m_csr_create(20,20,400,offsets,indices,values,&a)==M_OK);
 for(int repeat=0;repeat<100;repeat++){assert(m_csr_amd(&a,order,20)==M_OK);for(size_t i=0;i<20;i++)assert(order[i]==i);}
 m_csr_free(&a);return 0;
}
