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

#endif
