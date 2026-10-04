#ifndef LINEAR_A_MATRIX_H
#define LINEAR_A_MATRIX_H

#include <stdbool.h>
#include <stddef.h>

typedef struct { size_t rows, cols; double *values; } matrix;
typedef enum {
    M_OK, M_ARGUMENT, M_SHAPE, M_INDEX, M_MEMORY, M_NONFINITE,
    M_ALGORITHM, M_NOT_POSITIVE_DEFINITE, M_NOT_SYMMETRIC, M_NO_CONVERGENCE,
    M_SINGULAR, M_RANK_DEFICIENT, M_SOLVER_RANGE
} matrix_status;
typedef enum {
    M_DETERMINANT_AUTO = 0,
    M_DETERMINANT_LU = 1,
    M_DETERMINANT_CHOLESKY = 2
} matrix_determinant_algorithm;

typedef enum { M_LU=1, M_CHOLESKY=2, M_QR=3 } matrix_factor_algorithm;
typedef struct matrix_factor matrix_factor;
/* Factors own a snapshot; initialize *out to NULL, then free with m_factor_free.
 * QR uses column pivoting and requires rows >= cols and numerical full rank.
 * LU requires square input and nonzero computed pivots; Cholesky additionally
 * requires exact symmetry and positive computed pivots. Inputs are unchanged. */
matrix_status m_factorize(const matrix *source, matrix_factor_algorithm algorithm, matrix_factor **out);
void m_factor_free(matrix_factor **factor);
matrix_status m_factor_solve(const matrix_factor *factor, const matrix *rhs, matrix *out);
/* Computed infinity-norm reciprocal condition for square factors, O(n^3).
 * Diagnostic only, not a certified bound. Empty returns 1; inverse overflow 0. */
matrix_status m_factor_rcond(const matrix_factor *factor, double *out);
matrix_status m_solve(const matrix *source, const matrix *rhs, matrix *out);
matrix_status m_least_squares(const matrix *source, const matrix *rhs, matrix *out);

/* Initialize every output with {0}; free it before reusing it as an output.
 * Successful outputs own their allocation. Errors leave outputs unchanged.
 * Inputs and outputs must be distinct. Use m_set to update finite values. */
const char *m_error(matrix_status status);
matrix_status m_create(size_t rows, size_t cols, matrix *out);
matrix_status m_from_array(size_t rows, size_t cols, const double *values, matrix *out);
matrix_status m_identity(size_t size, matrix *out);
/* Right-handed active rotations of column vectors; angles are in radians.
 * Axis-angle accepts a nonzero finite 3-by-1 or 1-by-3 axis and normalizes it. */
matrix_status m_rotation_2d(double radians, matrix *out);
matrix_status m_rotation_x(double radians, matrix *out);
matrix_status m_rotation_y(double radians, matrix *out);
matrix_status m_rotation_z(double radians, matrix *out);
matrix_status m_rotation_axis_angle(const matrix *axis, double radians, matrix *out);
matrix_status m_copy(const matrix *source, matrix *out);
void m_free(matrix *value);
matrix_status m_get(const matrix *value, size_t row, size_t col, double *out);
matrix_status m_set(matrix *value, size_t row, size_t col, double number);
matrix_status m_add(const matrix *left, const matrix *right, matrix *out);
matrix_status m_subtract(const matrix *left, const matrix *right, matrix *out);
matrix_status m_scale(const matrix *source, double scalar, matrix *out);
matrix_status m_transpose(const matrix *source, matrix *out);
matrix_status m_multiply(const matrix *left, const matrix *right, matrix *out);
/* Three-dimensional vector cross product. Each operand may be 3-by-1 or
 * 1-by-3; the independently owned result retains the left operand's shape. */
matrix_status m_cross(const matrix *left, const matrix *right, matrix *out);
/* Extraction returns independently owned 1-by-cols / rows-by-1 matrices. */
matrix_status m_row(const matrix *source, size_t row, matrix *out);
matrix_status m_column(const matrix *source, size_t col, matrix *out);
matrix_status m_trace(const matrix *source, double *out);
/* Auto uses guarded formulas for sizes 0..3 and exact triangular structure,
 * falling back to partial-pivot LU. Cholesky requires exact symmetry and
 * positive computed pivots; Auto never assumes positive definiteness.
 * Determinant products are scaled to avoid intermediate range loss, while
 * elimination can still overflow. Final underflow follows double rounding.
 * All algorithms preserve the input and leave *out unchanged on error. */
matrix_status m_determinant(const matrix *source, double *out);
matrix_status m_determinant_with_algorithm(const matrix *source,
                                         matrix_determinant_algorithm algorithm,
                                         double *out);
matrix_status m_triangular(const matrix *source, bool *upper, bool *lower);

/* Real symmetric eigendecomposition A V = V diag(values). Values are sorted
 * ascending in an n-by-1 matrix; unit eigenvectors are columns of n-by-n V.
 * Exact symmetry is required. Cyclic Jacobi iterations use a relative
 * Frobenius off-diagonal tolerance (default 1e-12, at most 50 sweeps).
 * Options require 0 < tolerance < 1 and 1 <= max_sweeps <= 10000.
 * Both outputs must be empty, distinct from each other and the source.
 * Nonconvergence and nonfinite eigenvalues fail without changing outputs. */
matrix_status m_eigen_symmetric(const matrix *source, matrix *values, matrix *vectors);
matrix_status m_eigen_symmetric_with_options(const matrix *source, double tolerance,
                                            size_t max_sweeps, matrix *values, matrix *vectors);

/* General real-square eigendecomposition, including complex conjugate pairs.
 * A*(Vr+i*Vi) = (Vr+i*Vi)*diag(real+i*imag). Values are n-by-1, sorted by
 * (real, imag); complex vector columns have unit norm. Defective matrices may
 * have dependent columns. Hessenberg double-shift QR uses machine-precision
 * deflation and at most max_iterations steps between root deflations
 * (default 1000, allowed 1..100000). Four outputs must be empty and distinct.
 * Inputs are preserved and all outputs remain unchanged on any failure. */
matrix_status m_eigen_general(const matrix *source, matrix *values_real, matrix *values_imag,
                              matrix *vectors_real, matrix *vectors_imag);
matrix_status m_eigen_general_with_options(const matrix *source, size_t max_iterations,
                                          matrix *values_real, matrix *values_imag,
                                          matrix *vectors_real, matrix *vectors_imag);

/* Economy SVD A=U diag(values) Vt, k=min(rows,cols). U is rows-by-k,
 * values is k-by-1, Vt is k-by-cols. Descending nonnegative singular values.
 * Initialize distinct outputs with {0}; free each with m_free. Input unchanged.
 * Use tolerance=1e-12, max_sweeps=100 for defaults. See ports/SVD.md. */
matrix_status m_svd(const matrix *source, double tolerance, size_t max_sweeps,
                    matrix *u, matrix *values, matrix *vt);

/* Truncated SVD inverse and minimum-norm least squares. Cutoff in [0,1],
 * strict s/s_max > cutoff; -1 selects max(rows,cols)*DBL_EPSILON.
 * Empty distinct outputs, unchanged on failure. See ports/SVD.md. */
typedef struct { size_t rank; double reciprocal_condition, retained_reciprocal_condition; } matrix_spectral_diagnostics;
matrix_status m_pseudoinverse(const matrix *source,double cutoff,matrix *out);
matrix_status m_solve_ridge(const matrix *source,const matrix *rhs,double lambda,matrix *out);
matrix_status m_solve_minimum_norm(const matrix *source,const matrix *rhs,double cutoff,matrix *out);
matrix_status m_spectral_diagnostics(const matrix *source,double cutoff,matrix_spectral_diagnostics *out);

#include "sparse_core.h"
typedef la_csr sparse_matrix;
typedef la_cg_result matrix_cg_result;
/* Owned canonical CSR. See ports/SPARSE.md. All output handles start at {0}.
 * CG reason: 0 converged, 1 iteration_limit, 2 breakdown, 3 nonfinite. */
matrix_status m_csr_create(size_t rows,size_t cols,size_t nnz,const size_t *offsets,const size_t *indices,const double *values,sparse_matrix *out);
matrix_status m_csr_rcm(const sparse_matrix *a,size_t *out,size_t count);
matrix_status m_csr_amd(const sparse_matrix *a,size_t *out,size_t count);
matrix_status m_csr_permute(const sparse_matrix *a,const size_t *p,size_t count,sparse_matrix *out);
matrix_status m_permute_vector(const size_t *p,size_t count,const double *x,bool inverse,double *out);
matrix_status m_csr_from_dense(const matrix *a,sparse_matrix *out);
void m_csr_free(sparse_matrix *a);
matrix_status m_csr_matvec(const sparse_matrix *a,const matrix *x,matrix *out);
matrix_status m_csr_cg(const sparse_matrix *a,const double *b,size_t count,double rtol,double atol,size_t limit,bool jacobi,bool capture,matrix_cg_result *out);
void m_cg_free(matrix_cg_result *result);
#include "gmres_core.h"
typedef la_gmres_result matrix_gmres_result;
matrix_status m_csr_gmres(const sparse_matrix *a,const double *b,size_t count,size_t restart,double rtol,double atol,size_t limit,bool jacobi,bool capture,matrix_gmres_result *out);
void m_gmres_free(matrix_gmres_result *r);
typedef la_ilu0 matrix_ilu0;
matrix_status m_ilu0_create(const sparse_matrix *a,matrix_ilu0 *out);
void m_ilu0_free(matrix_ilu0 *f);
matrix_status m_ilu0_apply(const matrix_ilu0 *f,const matrix *b,matrix *out);
matrix_status m_csr_gmres_preconditioned(const sparse_matrix *a,const double *b,size_t count,size_t restart,double rtol,double atol,size_t limit,bool jacobi,bool capture,const matrix_ilu0 *preconditioner,matrix_gmres_result *out);
#include "cholesky_core.h"
typedef la_cholesky_symbolic matrix_cholesky_symbolic;
typedef la_cholesky matrix_cholesky;
matrix_status m_cholesky_analyze(const sparse_matrix *a,matrix_cholesky_symbolic *out);
void m_cholesky_symbolic_free(matrix_cholesky_symbolic *s);
matrix_status m_cholesky_factorize(const matrix_cholesky_symbolic *s,const sparse_matrix *a,matrix_cholesky *out);
void m_cholesky_free(matrix_cholesky *f);
matrix_status m_cholesky_solve(const matrix_cholesky *f,const matrix *b,matrix *out);
/* IC(0) retains the original lower pattern; nonpositive pivots fail without shifts. */
typedef la_ic0 matrix_ic0;
matrix_status m_ic0_create(const sparse_matrix *a,matrix_ic0 *out);
void m_ic0_free(matrix_ic0 *f);
matrix_status m_ic0_apply(const matrix_ic0 *f,const matrix *b,matrix *out);
matrix_status m_csr_cg_preconditioned(const sparse_matrix *a,const double *b,size_t count,double rtol,double atol,size_t limit,bool jacobi,bool capture,const matrix_ic0 *f,matrix_cg_result *out);

#include "multigrid_core.h"
typedef la_multigrid matrix_multigrid;
matrix_status m_multigrid_create(size_t width,matrix_multigrid *out);
matrix_status m_multigrid_matrix(const matrix_multigrid *m,sparse_matrix *out);
matrix_status m_multigrid_apply(const matrix_multigrid *m,const matrix *b,matrix *out);
matrix_status m_csr_cg_multigrid(const sparse_matrix *a,const double *b,size_t count,double rtol,double atol,size_t limit,bool jacobi,bool capture,const matrix_multigrid *m,matrix_cg_result *out);
#endif
