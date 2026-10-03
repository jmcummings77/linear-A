#include "matrix.h"
#include <assert.h>
#include <math.h>
int main(void){
 size_t rp[]={0,2,3},ci[]={0,1,1};double v[]={3,1,2},b[]={5,4};sparse_matrix a={0};matrix_gmres_result r={0};
 assert(m_csr_create(2,2,3,rp,ci,v,&a)==M_OK);
 assert(m_csr_gmres(&a,b,2,1,1e-10,0,100,0,1,&r)==M_OK);
 assert(r.reason==0&&r.restart_count>0&&fabs(r.x[0]-1)<1e-8&&fabs(r.x[1]-2)<1e-8);
 assert(r.residuals[r.iterations]<1e-9&&r.estimated_residuals[0]==r.residuals[0]);
 assert(m_csr_gmres(&a,b,2,1,1e-10,0,100,0,1,&r)==M_ARGUMENT);
 m_gmres_free(&r);m_gmres_free(&r);
 assert(m_csr_gmres(&a,b,2,0,1e-10,0,100,0,0,&r)==M_ARGUMENT);
 assert(m_csr_gmres(&a,b,2,2,1e-10,0,0,0,1,&r)==M_OK);assert(r.reason==1&&r.iterations==0&&r.iterates[0]==0);m_gmres_free(&r);
 m_csr_free(&a);
 size_t empty[]={0};assert(m_csr_create(0,0,0,empty,NULL,NULL,&a)==M_OK);
 assert(m_csr_gmres(&a,NULL,0,2,1e-10,0,10,0,1,&r)==M_OK);assert(r.reason==0);m_gmres_free(&r);m_csr_free(&a);
 return 0;
}
