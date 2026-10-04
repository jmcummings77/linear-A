#include "sparse_core.h"
static matrix_status sparse_status(int code){return code==0?M_OK:code==2?M_MEMORY:code==3?M_SOLVER_RANGE:code==4?M_NOT_SYMMETRIC:code==5?M_NOT_POSITIVE_DEFINITE:M_ARGUMENT;}
matrix_status m_csr_create(size_t rows,size_t cols,size_t nnz,const size_t *offsets,const size_t *indices,const double *values,sparse_matrix *out){return sparse_status(la_csr_create(rows,cols,nnz,offsets,indices,values,out));}
void m_csr_free(sparse_matrix *a){la_csr_free(a);}
void m_cg_free(matrix_cg_result *r){la_cg_free(r);}
matrix_status m_csr_from_dense(const matrix *a,sparse_matrix *out){
    matrix_status status=validate(a);if(status)return status;
    if(a->rows==SIZE_MAX||!la_sparse_count(a->rows+1,sizeof(size_t)))return M_MEMORY;
    size_t nnz=0;for(size_t i=0;i<a->rows*a->cols;i++)if(a->values[i]!=0)nnz++;
    size_t *rp=calloc(a->rows+1,sizeof(size_t)),*ci=calloc(nnz?nnz:1,sizeof(size_t));double *v=calloc(nnz?nnz:1,sizeof(double));
    if(!rp||!ci||!v){free(rp);free(ci);free(v);return M_MEMORY;}
    size_t k=0;for(size_t i=0;i<a->rows;i++){for(size_t j=0;j<a->cols;j++)if(a->values[i*a->cols+j]!=0){ci[k]=j;v[k++]=a->values[i*a->cols+j];}rp[i+1]=k;}
    status=m_csr_create(a->rows,a->cols,nnz,rp,ci,v,out);free(rp);free(ci);free(v);return status;
}
matrix_status m_csr_matvec(const sparse_matrix *a,const matrix *x,matrix *out){
    if(la_csr_validate(a)||!empty_output(out)||x==out)return M_ARGUMENT;
    matrix_status status=validate(x);if(status)return status;
    if(x->cols!=1||x->rows!=a->cols)return M_SHAPE;
    matrix result={0};status=m_create(a->rows,1,&result);if(status)return status;
    status=sparse_status(la_csr_mv(a,x->values,result.values));if(status)m_free(&result);else *out=result;return status;
}
matrix_status m_csr_cg(const sparse_matrix *a,const double *b,size_t count,double rtol,double atol,size_t limit,bool jacobi,bool capture,matrix_cg_result *out){return sparse_status(la_csr_cg(a,b,count,rtol,atol,limit,jacobi,capture,out));}

matrix_status m_csr_gmres(const sparse_matrix *a,const double *b,size_t count,size_t restart,double rtol,double atol,size_t limit,bool jacobi,bool capture,matrix_gmres_result *out){return sparse_status(la_csr_gmres(a,b,count,restart,rtol,atol,limit,jacobi,capture,out));}
void m_gmres_free(matrix_gmres_result *r){la_gmres_free(r);}

matrix_status m_ilu0_create(const sparse_matrix *a,matrix_ilu0 *out){int c=la_ilu0_create(a,out);return c==6?M_SINGULAR:sparse_status(c);}
void m_ilu0_free(matrix_ilu0 *f){la_ilu0_free(f);}
matrix_status m_ilu0_apply(const matrix_ilu0 *f,const matrix *b,matrix *out){
 if(la_ilu0_validate(f)||validate(b)||!empty_output(out)||b==out)return M_ARGUMENT;
 if(b->rows!=f->factors.rows||b->cols!=1)return M_SHAPE;
 matrix result={0};matrix_status status=m_create(b->rows,1,&result);if(status)return status;
 status=sparse_status(la_ilu0_apply(f,b->values,b->rows,result.values));if(status)m_free(&result);else *out=result;return status;
}
matrix_status m_csr_gmres_preconditioned(const sparse_matrix *a,const double *b,size_t count,size_t restart,double rtol,double atol,size_t limit,bool jacobi,bool capture,const matrix_ilu0 *preconditioner,matrix_gmres_result *out){return sparse_status(la_csr_gmres_preconditioned(a,b,count,restart,rtol,atol,limit,jacobi,capture,preconditioner,out));}

#include "ordering_core.h"
matrix_status m_csr_rcm(const sparse_matrix *a,size_t *out,size_t count){return sparse_status(la_csr_rcm(a,out,count));}
matrix_status m_csr_amd(const sparse_matrix *a,size_t *out,size_t count){return sparse_status(la_csr_amd(a,out,count));}
matrix_status m_csr_permute(const sparse_matrix *a,const size_t *p,size_t count,sparse_matrix *out){return sparse_status(la_csr_permute(a,p,count,out));}
matrix_status m_permute_vector(const size_t *p,size_t count,const double *x,bool inverse,double *out){return sparse_status(la_permute_vector(p,count,x,inverse,out));}

matrix_status m_cholesky_analyze(const sparse_matrix *a,matrix_cholesky_symbolic *out){return sparse_status(la_cholesky_analyze(a,out));}
void m_cholesky_symbolic_free(matrix_cholesky_symbolic *s){la_cholesky_symbolic_free(s);}
matrix_status m_cholesky_factorize(const matrix_cholesky_symbolic *s,const sparse_matrix *a,matrix_cholesky *out){return sparse_status(la_cholesky_factorize(s,a,out));}
void m_cholesky_free(matrix_cholesky *f){la_cholesky_free(f);}
matrix_status m_cholesky_solve(const matrix_cholesky *f,const matrix *b,matrix *out){
 if(!f||validate(b)||!empty_output(out)||b==out)return M_ARGUMENT;
 if(b->cols!=1||b->rows!=f->lower.rows)return M_SHAPE;
 matrix result={0};matrix_status code=m_create(b->rows,1,&result);if(code)return code;
 code=sparse_status(la_cholesky_solve(f,b->values,b->rows,result.values));if(code)m_free(&result);else *out=result;return code;
}

matrix_status m_ic0_create(const sparse_matrix *a,matrix_ic0 *out){return sparse_status(la_ic0_create(a,out));}
void m_ic0_free(matrix_ic0 *f){la_cholesky_free(f);}
matrix_status m_ic0_apply(const matrix_ic0 *f,const matrix *b,matrix *out){return m_cholesky_solve(f,b,out);}
matrix_status m_csr_cg_preconditioned(const sparse_matrix *a,const double *b,size_t count,double rtol,double atol,size_t limit,bool jacobi,bool capture,const matrix_ic0 *f,matrix_cg_result *out){return sparse_status(la_csr_cg_ic0(a,b,count,rtol,atol,limit,jacobi,capture,f,out));}

matrix_status m_multigrid_create(size_t width,matrix_multigrid *out){return sparse_status(la_mg_create(width,out));}
matrix_status m_multigrid_matrix(const matrix_multigrid *m,sparse_matrix *out){return sparse_status(la_mg_matrix(m,out));}
matrix_status m_multigrid_apply(const matrix_multigrid *m,const matrix *b,matrix *out){if(!la_mg_valid(m)||validate(b)||!empty_output(out)||b==out)return M_ARGUMENT;if(b->cols!=1||b->rows!=m->width*m->width)return M_SHAPE;matrix result={0};matrix_status code=m_create(b->rows,1,&result);if(code)return code;code=sparse_status(la_mg_apply(m,b->values,b->rows,result.values));if(code)m_free(&result);else *out=result;return code;}
matrix_status m_csr_cg_multigrid(const sparse_matrix *a,const double *b,size_t count,double rtol,double atol,size_t limit,bool jacobi,bool capture,const matrix_multigrid *m,matrix_cg_result *out){return sparse_status(la_csr_cg_mg(a,b,count,rtol,atol,limit,jacobi,capture,m,out));}
