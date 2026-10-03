/* Opaque owned handles around the reusable C matrix implementation. */
#include "../c/matrix.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#ifdef MATRIX_TRACE
#include "trace.h"
#endif

static matrix_status last_status = M_OK;

int wm_last_error(void) { return last_status; }
const char *wm_error_message(void) { return m_error(last_status); }

static matrix *allocate_handle(void) {
    matrix *result = calloc(1, sizeof(matrix));
    if (!result) last_status = M_MEMORY;
    return result;
}

static matrix *complete(matrix *result, matrix_status status) {
    last_status = status;
    if (status != M_OK) { m_free(result); free(result); return NULL; }
    return result;
}

matrix *wm_create(size_t rows, size_t cols) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_create(rows, cols, result)) : NULL;
}

matrix *wm_identity(size_t size) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_identity(size, result)) : NULL;
}

void wm_destroy(matrix *value) { if (value) { m_free(value); free(value); } }
size_t wm_rows(const matrix *value) { return value->rows; }
size_t wm_cols(const matrix *value) { return value->cols; }
double *wm_data(matrix *value) { return value->values; }

matrix *wm_copy(const matrix *source) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_copy(source, result)) : NULL;
}

matrix *wm_add(const matrix *a, const matrix *b) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_add(a, b, result)) : NULL;
}

matrix *wm_subtract(const matrix *a, const matrix *b) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_subtract(a, b, result)) : NULL;
}

matrix *wm_scale(const matrix *a, double scalar) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_scale(a, scalar, result)) : NULL;
}

matrix *wm_transpose(const matrix *a) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_transpose(a, result)) : NULL;
}

matrix *wm_multiply(const matrix *a, const matrix *b) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_multiply(a, b, result)) : NULL;
}

matrix *wm_cross(const matrix *a, const matrix *b) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_cross(a, b, result)) : NULL;
}

matrix *wm_rotation_2d(double radians) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_rotation_2d(radians, result)) : NULL;
}

matrix *wm_rotation_x(double radians) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_rotation_x(radians, result)) : NULL;
}

matrix *wm_rotation_y(double radians) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_rotation_y(radians, result)) : NULL;
}

matrix *wm_rotation_z(double radians) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_rotation_z(radians, result)) : NULL;
}

matrix *wm_rotation_axis_angle(const matrix *axis, double radians) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_rotation_axis_angle(axis, radians, result)) : NULL;
}

#ifdef MATRIX_TRACE
matrix *wm_multiply_trace(const matrix *a, const matrix *b) {
    m_trace_clear();
    if (!a || !b) { last_status = M_ARGUMENT; return NULL; }
    if (a->cols != b->rows) { last_status = M_SHAPE; return NULL; }
    if (a->rows > M_TRACE_MAX_DIMENSION || a->cols > M_TRACE_MAX_DIMENSION ||
        b->rows > M_TRACE_MAX_DIMENSION || b->cols > M_TRACE_MAX_DIMENSION) {
        last_status = M_ARGUMENT;
        return NULL;
    }
    matrix *result = allocate_handle();
    if (!result) return NULL;
    m_trace_begin();
    matrix_status status = m_multiply(a, b, result);
    if (!m_trace_end(status == M_OK)) status = M_MEMORY;
    return complete(result, status);
}

size_t wm_trace_count(void) { return m_trace_count(); }
const m_trace_step *wm_trace_data(void) { return m_trace_data(); }
size_t wm_trace_stride(void) { return sizeof(m_trace_step) / sizeof(double); }
#endif

matrix *wm_row(const matrix *a, size_t index) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_row(a, index, result)) : NULL;
}

matrix *wm_column(const matrix *a, size_t index) {
    matrix *result = allocate_handle();
    return result ? complete(result, m_column(a, index, result)) : NULL;
}

int wm_set(matrix *a, size_t row, size_t col, double value) {
    return last_status = m_set(a, row, col, value);
}

double wm_trace(const matrix *a) {
    double value = NAN;
    last_status = m_trace(a, &value);
    return value;
}

double wm_determinant(const matrix *a) {
    double value = NAN;
    last_status = m_determinant(a, &value);
    return value;
}

double wm_determinant_algorithm(const matrix *a, int algorithm) {
    double value = NAN;
    last_status = m_determinant_with_algorithm(a, (matrix_determinant_algorithm)algorithm, &value);
    return value;
}

/* A temporary packed handle gives JS one allocation to own across the call.
 * Row i contains eigenvalue i, followed by row i of the eigenvector matrix. */
matrix *wm_eigen_symmetric(const matrix *a, double tolerance, size_t max_sweeps) {
    matrix values = {0}, vectors = {0};
    last_status = m_eigen_symmetric_with_options(a, tolerance, max_sweeps, &values, &vectors);
    if (last_status != M_OK) return NULL;
    matrix *result = allocate_handle();
    if (!result) { m_free(&values); m_free(&vectors); return NULL; }
    const size_t n = values.rows;
    last_status = n == SIZE_MAX ? M_MEMORY : m_create(n, n + 1, result);
    if (last_status == M_OK) for (size_t row = 0; row < n; row++) {
        result->values[row * (n + 1)] = values.values[row];
        for (size_t col = 0; col < n; col++) result->values[row * (n + 1) + col + 1] = vectors.values[row * n + col];
    }
    m_free(&values); m_free(&vectors);
    return complete(result, last_status);
}

/* Row i packs real/imag values, followed by real/imag eigenvector rows. */
matrix *wm_eigen_general(const matrix *a, size_t max_iterations) {
    matrix dr = {0}, di = {0}, vr = {0}, vi = {0};
    last_status = m_eigen_general_with_options(a, max_iterations, &dr, &di, &vr, &vi);
    if (last_status != M_OK) return NULL;
    matrix *result = allocate_handle();
    if (!result) { m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi); return NULL; }
    const size_t n = dr.rows;
    last_status = n > (SIZE_MAX - 2) / 2 ? M_MEMORY : m_create(n, 2 * n + 2, result);
    if (last_status == M_OK) for (size_t row = 0; row < n; row++) {
        const size_t offset = row * (2 * n + 2);
        result->values[offset] = dr.values[row]; result->values[offset + 1] = di.values[row];
        for (size_t col = 0; col < n; col++) {
            result->values[offset + 2 + col] = vr.values[row * n + col];
            result->values[offset + 2 + n + col] = vi.values[row * n + col];
        }
    }
    m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi);
    return complete(result, last_status);
}

int wm_triangular(const matrix *a) {
    bool upper = false, lower = false;
    last_status = m_triangular(a, &upper, &lower);
    return (upper ? 1 : 0) | (lower ? 2 : 0);
}

double wm_checksum(const matrix *a) {
    const size_t count = a->rows * a->cols;
    double result = count ? a->values[0] + a->values[count / 2] + a->values[count - 1] : 0.0;
    last_status = isfinite(result) ? M_OK : M_NONFINITE;
    return result;
}

matrix_factor *wm_factorize(const matrix *source,int algorithm){
    matrix_factor *factor=NULL;last_status=m_factorize(source,(matrix_factor_algorithm)algorithm,&factor);return factor;
}
void wm_factor_destroy(matrix_factor *factor){m_factor_free(&factor);}
matrix *wm_factor_solve(const matrix_factor *factor,const matrix *rhs){
    matrix *result=allocate_handle();return result?complete(result,m_factor_solve(factor,rhs,result)):NULL;
}
double wm_factor_rcond(const matrix_factor *factor){double result=0;last_status=m_factor_rcond(factor,&result);return result;}

int wm_svd(const matrix *source,double tolerance,size_t max_sweeps,matrix *u,matrix *values,matrix *vt){
    last_status=m_svd(source,tolerance,max_sweeps,u,values,vt);
    return (int)last_status;
}

matrix *wm_pseudoinverse(const matrix *source,double cutoff){
    matrix *out=allocate_handle();return out?complete(out,m_pseudoinverse(source,cutoff,out)):NULL;
}
matrix *wm_solve_minimum_norm(const matrix *source,const matrix *rhs,double cutoff){
    matrix *out=allocate_handle();return out?complete(out,m_solve_minimum_norm(source,rhs,cutoff,out)):NULL;
}
matrix *wm_spectral_diagnostics(const matrix *source,double cutoff){
    matrix_spectral_diagnostics d;
    last_status=m_spectral_diagnostics(source,cutoff,&d);if(last_status!=M_OK)return NULL;
    matrix *out=allocate_handle();if(!out)return NULL;
    double values[3]={(double)d.rank,d.reciprocal_condition,d.retained_reciprocal_condition};
    return complete(out,m_from_array(1,3,values,out));
}

matrix *wm_solve_ridge(const matrix *a,const matrix *b,double lambda){
    matrix *out=allocate_handle();return out?complete(out,m_solve_ridge(a,b,lambda,out)):NULL;
}

sparse_matrix *wm_csr_create(size_t rows,size_t cols,const matrix *rp,const matrix *ci,const matrix *v){
    last_status=M_ARGUMENT;
    if(!rp||!ci||!v||rows==SIZE_MAX||rp->rows!=rows+1||rp->cols!=1||ci->cols!=1||v->cols!=1||ci->rows!=v->rows)return NULL;
    size_t *offsets=calloc(rows+1,sizeof(size_t)),*indices=calloc(ci->rows?ci->rows:1,sizeof(size_t));
    sparse_matrix *out=calloc(1,sizeof(sparse_matrix));
    if(!offsets||!indices||!out){free(offsets);free(indices);free(out);last_status=M_MEMORY;return NULL;}
    for(size_t i=0;i<=rows;i++){double x=rp->values[i];if(!isfinite(x)||x<0||x>SIZE_MAX||floor(x)!=x)goto invalid;offsets[i]=(size_t)x;}
    for(size_t i=0;i<ci->rows;i++){double x=ci->values[i];if(!isfinite(x)||x<0||x>SIZE_MAX||floor(x)!=x)goto invalid;indices[i]=(size_t)x;}
    last_status=m_csr_create(rows,cols,v->rows,offsets,indices,v->values,out);free(offsets);free(indices);
    if(last_status){free(out);return NULL;}return out;
invalid:free(offsets);free(indices);free(out);return NULL;
}
void wm_csr_destroy(sparse_matrix *a){if(a){m_csr_free(a);free(a);}}
matrix *wm_csr_matvec(const sparse_matrix *a,const matrix *x){matrix *out=allocate_handle();return out?complete(out,m_csr_matvec(a,x,out)):NULL;}
matrix_cg_result *wm_csr_cg(const sparse_matrix *a,const matrix *b,double rtol,double atol,size_t limit,int jacobi,int capture){
    if(!b||b->cols!=1){last_status=M_ARGUMENT;return NULL;}
    matrix_cg_result *out=calloc(1,sizeof(matrix_cg_result));if(!out){last_status=M_MEMORY;return NULL;}
    last_status=m_csr_cg(a,b->values,b->rows,rtol,atol,limit,jacobi!=0,capture!=0,out);
    if(last_status){free(out);return NULL;}return out;
}
void wm_cg_destroy(matrix_cg_result *r){if(r){m_cg_free(r);free(r);}}
size_t wm_cg_iterations(const matrix_cg_result *r){return r->iterations;}
int wm_cg_reason(const matrix_cg_result *r){return r->reason;}
double *wm_cg_data(matrix_cg_result *r,int field){return field==0?r->x:field==1?r->residuals:r->iterates;}

matrix_gmres_result *wm_csr_gmres(const sparse_matrix *a,const matrix *b,size_t restart,double rtol,double atol,size_t limit,int jacobi,int capture){
    if(!b||b->cols!=1){last_status=M_ARGUMENT;return NULL;}
    matrix_gmres_result *out=calloc(1,sizeof(matrix_gmres_result));if(!out){last_status=M_MEMORY;return NULL;}
    last_status=m_csr_gmres(a,b->values,b->rows,restart,rtol,atol,limit,jacobi!=0,capture!=0,out);
    if(last_status){free(out);return NULL;}return out;
}
void wm_gmres_destroy(matrix_gmres_result *r){if(r){m_gmres_free(r);free(r);}}
size_t wm_gmres_iterations(const matrix_gmres_result *r){return r->iterations;}
int wm_gmres_reason(const matrix_gmres_result *r){return r->reason;}
size_t wm_gmres_restart_count(const matrix_gmres_result *r){return r->restart_count;}
size_t *wm_gmres_restarts(matrix_gmres_result *r){return r->restarts;}
double *wm_gmres_data(matrix_gmres_result *r,int field){return field==0?r->x:field==1?r->residuals:field==2?r->estimated_residuals:r->iterates;}
