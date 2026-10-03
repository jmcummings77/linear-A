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
