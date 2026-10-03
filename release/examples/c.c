#include "matrix.h"
#include <math.h>
#include <stdio.h>
int main(void) {
    matrix a={0}, b={0}, x={0}, reconstructed={0};
    double av[]={4,1,2,3}, bv[]={6,8};
    int failed=m_from_array(2,2,av,&a) || m_from_array(2,1,bv,&b)
        || m_solve(&a,&b,&x) || m_multiply(&a,&x,&reconstructed);
    if (!failed) for (size_t i=0;i<2;i++)
        failed |= !isfinite(x.values[i]) || fabs(x.values[i]-(i+1))>1e-12
            || !isfinite(reconstructed.values[i]) || fabs(reconstructed.values[i]-bv[i])>1e-12;
    matrix u={0},s={0},vt={0};
    if (!failed) failed=m_svd(&a,1e-12,100,&u,&s,&vt);
    if (!failed) failed=s.rows!=2 || s.values[1]<=0;
    matrix inverse={0},minimum={0};matrix_spectral_diagnostics diagnostics;
    if(!failed)failed=m_pseudoinverse(&a,-1,&inverse)||m_solve_minimum_norm(&a,&b,-1,&minimum)||m_spectral_diagnostics(&a,-1,&diagnostics);
    if(!failed)failed=fabs(inverse.values[0]-.3)>1e-12||fabs(minimum.values[1]-2)>1e-12||diagnostics.rank!=2;
    m_free(&inverse);m_free(&minimum);
    m_free(&u);m_free(&s);m_free(&vt);
    matrix ridge={0};
    if(!failed)failed=m_solve_ridge(&a,&b,1,&ridge);
    if(!failed)failed=fabs(ridge.values[0]-140.0/131)>1e-12||fabs(ridge.values[1]-230.0/131)>1e-12;
    m_free(&ridge);
    m_free(&reconstructed);m_free(&x);m_free(&b);m_free(&a);
    sparse_matrix sparse={0};matrix_cg_result cg={0};size_t rp[]={0,2,4},ci[]={0,1,0,1};double sv[]={4,1,1,3},rhs[]={6,7};
    if(!failed)failed=m_csr_create(2,2,4,rp,ci,sv,&sparse)||m_csr_cg(&sparse,rhs,2,1e-10,0,1000,true,true,&cg);
    if(!failed)failed=cg.reason!=0||fabs(cg.x[0]-1)>1e-12||fabs(cg.x[1]-2)>1e-12;
    matrix vector={0},product={0};double vector_values[]={1,2};
    if(!failed)failed=m_from_array(2,1,vector_values,&vector)||m_csr_matvec(&sparse,&vector,&product);
    if(!failed)failed=product.values[0]!=6||product.values[1]!=7;
    matrix_gmres_result gm={0};
    if(!failed)failed=m_csr_gmres(&sparse,rhs,2,2,1e-10,0,1000,true,true,&gm);
    if(!failed)failed=gm.reason!=0||fabs(gm.x[0]-1)>1e-12||fabs(gm.x[1]-2)>1e-12;
    m_gmres_free(&gm);
    matrix_ilu0 ilu={0};matrix applied={0};
    if(!failed)failed=m_ilu0_create(&sparse,&ilu);
    if(!failed){vector.values[0]=6;vector.values[1]=7;failed=m_ilu0_apply(&ilu,&vector,&applied);}
    if(!failed)failed=fabs(applied.values[0]-1)>1e-12;m_free(&applied);
    if(!failed){vector.values[0]=11;vector.values[1]=13;failed=m_ilu0_apply(&ilu,&vector,&applied);}
    if(!failed)failed=fabs(applied.values[0]-20.0/11)>1e-12;m_free(&applied);
    if(!failed)failed=m_csr_gmres_preconditioned(&sparse,rhs,2,20,1e-10,0,1000,false,false,&ilu,&gm);
    if(!failed)failed=gm.reason!=0;m_gmres_free(&gm);m_ilu0_free(&ilu);
    size_t order[2];sparse_matrix reordered={0};double y[2],restored[2],input[]={1,2};
    if(!failed)failed=m_csr_rcm(&sparse,order,2)||m_csr_permute(&sparse,order,2,&reordered)||m_permute_vector(order,2,input,false,y)||m_permute_vector(order,2,y,true,restored);
    if(!failed)failed=restored[0]!=1||restored[1]!=2;m_csr_free(&reordered);
    matrix_cholesky_symbolic plan={0};matrix_cholesky chol={0};matrix cx={0};
    if(!failed)failed=m_cholesky_analyze(&sparse,&plan)||m_cholesky_factorize(&plan,&sparse,&chol);
    if(!failed){vector.values[0]=6;vector.values[1]=7;failed=m_cholesky_solve(&chol,&vector,&cx);}
    if(!failed)failed=fabs(cx.values[0]-1)>1e-12;m_free(&cx);
    if(!failed){vector.values[0]=11;vector.values[1]=13;failed=m_cholesky_solve(&chol,&vector,&cx);}
    if(!failed)failed=fabs(cx.values[0]-20.0/11)>1e-12;m_free(&cx);m_cholesky_free(&chol);m_cholesky_symbolic_free(&plan);
    m_free(&vector);m_free(&product);m_cg_free(&cg);m_csr_free(&sparse);
    if (!failed) puts("solution: 1, 2");
    return failed;
}
