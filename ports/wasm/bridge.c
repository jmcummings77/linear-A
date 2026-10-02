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
