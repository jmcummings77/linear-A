#include "matrix.h"
#include <float.h>
#include <limits.h>
#include "general_eigen.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef MATRIX_TRACE
#include "../wasm/trace.h"
#ifdef MATRIX_USE_ASM
#error MATRIX_TRACE captures the C multiplication loop and cannot use assembly kernels.
#endif
#endif

#ifdef MATRIX_USE_ASM
#if !defined(__aarch64__) && !defined(__arm64__)
#error MATRIX_USE_ASM requires ARM64; there is no C fallback.
#endif
extern void m_asm_add(const double *, const double *, double *, size_t);
extern void m_asm_subtract(const double *, const double *, double *, size_t);
extern void m_asm_scale(const double *, double *, size_t, double);
extern void m_asm_transpose(const double *, double *, size_t, size_t);
extern void m_asm_multiply(const double *, const double *, double *, size_t, size_t, size_t);
extern double m_asm_trace(const double *, size_t);
extern void m_asm_eliminate(double *, const double *, size_t, double);
extern void m_asm_cross(const double *, const double *, double *);
#endif

const char *m_error(matrix_status status) {
    switch (status) {
        case M_OK: return "success";
        case M_ARGUMENT: return "invalid argument or nonempty/aliased output";
        case M_SHAPE: return "incompatible matrix dimensions";
        case M_INDEX: return "matrix index out of bounds";
        case M_MEMORY: return "allocation failed or dimensions overflow";
        case M_NONFINITE: return "nonfinite input or arithmetic result";
        case M_ALGORITHM: return "unknown determinant algorithm";
        case M_NOT_POSITIVE_DEFINITE: return "Cholesky requires an exactly symmetric positive definite matrix";
        case M_NOT_SYMMETRIC: return "eigendecomposition requires an exactly symmetric matrix";
        case M_NO_CONVERGENCE: return "eigendecomposition did not converge";
        case M_SINGULAR: return "singular matrix: zero computed LU pivot";
        case M_RANK_DEFICIENT: return "QR input is numerically rank deficient";
        case M_SOLVER_RANGE: return "solver scaling or arithmetic exceeds the float64 range";
    }
    return "unknown matrix error";
}

static bool empty_output(const matrix *out) {
    return out && !out->values && !out->rows && !out->cols;
}

static matrix_status validate(const matrix *m) {
    if (!m) return M_ARGUMENT;
    if (m->cols && m->rows > SIZE_MAX / m->cols / sizeof(double)) return M_MEMORY;
    size_t count = m->rows * m->cols;
    if (count && !m->values) return M_ARGUMENT;
    for (size_t i = 0; i < count; i++) if (!isfinite(m->values[i])) return M_NONFINITE;
    return M_OK;
}

matrix_status m_create(size_t rows, size_t cols, matrix *out) {
    if (!empty_output(out)) return M_ARGUMENT;
    if (cols && rows > SIZE_MAX / cols / sizeof(double)) return M_MEMORY;
    size_t count = rows * cols;
    double *values = count ? calloc(count, sizeof(double)) : NULL;
    if (count && !values) return M_MEMORY;
    *out = (matrix){rows, cols, values};
    return M_OK;
}

void m_free(matrix *value) {
    if (value) { free(value->values); *value = (matrix){0}; }
}

matrix_status m_from_array(size_t rows, size_t cols, const double *values, matrix *out) {
    matrix result = {0};
    if (!empty_output(out) || (rows && cols && !values)) return M_ARGUMENT;
    matrix_status status = m_create(rows, cols, &result);
    if (status != M_OK) return status;
    for (size_t i = 0; i < rows * cols; i++) {
        if (!isfinite(values[i])) { m_free(&result); return M_NONFINITE; }
        result.values[i] = values[i];
    }
    *out = result;
    return M_OK;
}

matrix_status m_copy(const matrix *source, matrix *out) {
    if (source == out) return M_ARGUMENT;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    return m_from_array(source->rows, source->cols, source->values, out);
}

matrix_status m_identity(size_t size, matrix *out) {
    matrix_status status = m_create(size, size, out);
    if (status != M_OK) return status;
    for (size_t i = 0; i < size; i++) out->values[i * size + i] = 1.0;
    return M_OK;
}

matrix_status m_get(const matrix *value, size_t row, size_t col, double *out) {
    if (!value || !out) return M_ARGUMENT;
    if (row >= value->rows || col >= value->cols) return M_INDEX;
    *out = value->values[row * value->cols + col];
    return M_OK;
}

matrix_status m_set(matrix *value, size_t row, size_t col, double number) {
    if (!value) return M_ARGUMENT;
    if (row >= value->rows || col >= value->cols) return M_INDEX;
    if (!isfinite(number)) return M_NONFINITE;
    value->values[row * value->cols + col] = number;
    return M_OK;
}

static matrix_status finish(matrix *result, matrix *out) {
    matrix_status status = validate(result);
    if (status != M_OK) { m_free(result); return status; }
    *out = *result;
    return M_OK;
}

static matrix_status combine(const matrix *a, const matrix *b, matrix *out, bool subtract) {
    if (!empty_output(out) || a == out || b == out) return M_ARGUMENT;
    matrix_status status = validate(a);
    if (status != M_OK) return status;
    status = validate(b);
    if (status != M_OK) return status;
    if (a->rows != b->rows || a->cols != b->cols) return M_SHAPE;
    matrix result = {0};
    status = m_create(a->rows, a->cols, &result);
    if (status != M_OK) return status;
    size_t count = a->rows * a->cols;
#ifdef MATRIX_USE_ASM
    if (subtract) m_asm_subtract(a->values, b->values, result.values, count);
    else m_asm_add(a->values, b->values, result.values, count);
#else
    for (size_t i = 0; i < count; i++)
        result.values[i] = subtract ? a->values[i] - b->values[i] : a->values[i] + b->values[i];
#endif
    return finish(&result, out);
}

matrix_status m_add(const matrix *a, const matrix *b, matrix *out) { return combine(a, b, out, false); }
matrix_status m_subtract(const matrix *a, const matrix *b, matrix *out) { return combine(a, b, out, true); }

matrix_status m_scale(const matrix *source, double scalar, matrix *out) {
    if (!empty_output(out) || source == out) return M_ARGUMENT;
    if (!isfinite(scalar)) return M_NONFINITE;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    matrix result = {0};
    status = m_create(source->rows, source->cols, &result);
    if (status != M_OK) return status;
    size_t count = source->rows * source->cols;
#ifdef MATRIX_USE_ASM
    m_asm_scale(source->values, result.values, count, scalar);
#else
    for (size_t i = 0; i < count; i++) result.values[i] = source->values[i] * scalar;
#endif
    return finish(&result, out);
}

matrix_status m_transpose(const matrix *source, matrix *out) {
    if (!empty_output(out) || source == out) return M_ARGUMENT;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    matrix result = {0};
    status = m_create(source->cols, source->rows, &result);
    if (status != M_OK) return status;
#ifdef MATRIX_USE_ASM
    m_asm_transpose(source->values, result.values, source->rows, source->cols);
#else
    for (size_t row = 0; row < source->rows; row++)
    for (size_t col = 0; col < source->cols; col++)
        result.values[col * source->rows + row] = source->values[row * source->cols + col];
#endif
    return finish(&result, out);
}

matrix_status m_multiply(const matrix *a, const matrix *b, matrix *out) {
    if (!empty_output(out) || a == out || b == out) return M_ARGUMENT;
    matrix_status status = validate(a);
    if (status != M_OK) return status;
    status = validate(b);
    if (status != M_OK) return status;
    if (a->cols != b->rows) return M_SHAPE;
    matrix result = {0};
    status = m_create(a->rows, b->cols, &result);
    if (status != M_OK) return status;
#ifdef MATRIX_USE_ASM
    m_asm_multiply(a->values, b->values, result.values, a->rows, a->cols, b->cols);
#else
    for (size_t row = 0; row < a->rows; row++)
    for (size_t col = 0; col < b->cols; col++) {
        double sum = 0.0;
#ifdef MATRIX_TRACE
        for (size_t k = 0; k < a->cols; k++) {
            double left = a->values[row * a->cols + k];
            double right = b->values[k * b->cols + col];
            double product = left * right;
            sum += product;
            m_trace_record(row, col, k, left, right, product, sum, __LINE__);
        }
#else
        for (size_t k = 0; k < a->cols; k++) sum += a->values[row * a->cols + k] * b->values[k * b->cols + col];
#endif
        result.values[row * b->cols + col] = sum;
    }
#endif
    return finish(&result, out);
}

matrix_status m_row(const matrix *source, size_t row, matrix *out) {
    if (!empty_output(out) || source == out) return M_ARGUMENT;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    if (row >= source->rows) return M_INDEX;
    return m_from_array(1, source->cols, source->cols ? source->values + row * source->cols : NULL, out);
}

matrix_status m_column(const matrix *source, size_t col, matrix *out) {
    if (!empty_output(out) || source == out) return M_ARGUMENT;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    if (col >= source->cols) return M_INDEX;
    status = m_create(source->rows, 1, out);
    if (status != M_OK) return status;
    for (size_t row = 0; row < source->rows; row++) out->values[row] = source->values[row * source->cols + col];
    return M_OK;
}

matrix_status m_trace(const matrix *source, double *out) {
    if (!out) return M_ARGUMENT;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    if (source->rows != source->cols) return M_SHAPE;
    double sum = 0.0;
#ifdef MATRIX_USE_ASM
    sum = m_asm_trace(source->values, source->rows);
#else
    for (size_t i = 0; i < source->rows; i++) sum += source->values[i * source->cols + i];
#endif
    if (!isfinite(sum)) return M_NONFINITE;
    *out = sum;
    return M_OK;
}

/* Keep determinant factors normalized until the final result. A sequence of
 * large and small pivots must not overflow or vanish merely because of order. */
typedef struct { double fraction; int64_t exponent; } determinant_product;

static void product_append(determinant_product *product, double factor) {
    int exponent = 0, normalization = 0;
    double fraction = frexp(factor, &exponent);
    product->fraction = frexp(product->fraction * fraction, &normalization);
    product->exponent += exponent + normalization;
}

static matrix_status product_finish(determinant_product product, double *out) {
    double result;
    if (product.fraction == 0.0) result = product.fraction;
    else if (product.exponent > DBL_MAX_EXP) return M_NONFINITE;
    else if (product.exponent < DBL_MIN_EXP - DBL_MANT_DIG) result = copysign(0.0, product.fraction);
    else result = scalbn(product.fraction, (int)product.exponent);
    if (!isfinite(result)) return M_NONFINITE;
    *out = result;
    return M_OK;
}

static double quotient_product(double numerator, double denominator, double value) {
    int numerator_exponent, denominator_exponent, value_exponent;
    double fraction = frexp(numerator, &numerator_exponent) /
                      frexp(denominator, &denominator_exponent) *
                      frexp(value, &value_exponent);
    return scalbn(fraction, numerator_exponent - denominator_exponent + value_exponent);
}

static matrix_status determinant_lu(const matrix *source, double *out) {
    matrix work = {0};
    matrix_status status = m_copy(source, &work);
    if (status != M_OK) return status;
    size_t n = source->rows;
    determinant_product product = {1.0, 0};
    for (size_t k = 0; k < n; k++) {
        size_t pivot = k;
        for (size_t row = k + 1; row < n; row++)
            if (fabs(work.values[row * n + k]) > fabs(work.values[pivot * n + k])) pivot = row;
        if (work.values[pivot * n + k] == 0.0) { product.fraction = 0.0; break; }
        if (pivot != k) {
            for (size_t col = 0; col < n; col++) {
                double tmp = work.values[k * n + col];
                work.values[k * n + col] = work.values[pivot * n + col];
                work.values[pivot * n + col] = tmp;
            }
            product.fraction = -product.fraction;
        }
        double diagonal = work.values[k * n + k];
        product_append(&product, diagonal);
        for (size_t row = k + 1; row < n; row++) {
            double numerator = work.values[row * n + k];
            double factor = numerator / diagonal;
            work.values[row * n + k] = 0.0;
            /* A tiny ratio can lose range before multiplication by a large
             * pivot-row entry restores it. Preserve that update explicitly. */
            if (numerator != 0.0 && fabs(factor) < DBL_MIN) {
                for (size_t col = k + 1; col < n; col++)
                    work.values[row * n + col] -= quotient_product(numerator, diagonal, work.values[k * n + col]);
            } else {
#ifdef MATRIX_USE_ASM
                m_asm_eliminate(work.values + row * n + k + 1, work.values + k * n + k + 1, n - k - 1, factor);
#else
                for (size_t col = k + 1; col < n; col++)
                    work.values[row * n + col] -= factor * work.values[k * n + col];
#endif
            }
            for (size_t col = k + 1; col < n; col++)
                if (!isfinite(work.values[row * n + col])) { m_free(&work); return M_NONFINITE; }
        }
    }
    m_free(&work);
    return product_finish(product, out);
}

static matrix_status determinant_cholesky(const matrix *source, double *out) {
    size_t n = source->rows;
    for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < row; col++)
        if (source->values[row * n + col] != source->values[col * n + row])
            return M_NOT_POSITIVE_DEFINITE;
    matrix lower = {0};
    matrix_status status = m_create(n, n, &lower);
    if (status != M_OK) return status;
    determinant_product product = {1.0, 0};
    for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col <= row; col++) {
        double value = source->values[row * n + col];
        for (size_t k = 0; k < col; k++)
            value -= lower.values[row * n + k] * lower.values[col * n + k];
        if (!isfinite(value)) { m_free(&lower); return M_NONFINITE; }
        if (row == col) {
            if (value <= 0.0) { m_free(&lower); return M_NOT_POSITIVE_DEFINITE; }
            lower.values[row * n + col] = sqrt(value);
            /* The diagonal pivot is L[row,row]^2 before sqrt rounding. */
            product_append(&product, value);
        } else {
            value /= lower.values[col * n + col];
            if (!isfinite(value)) { m_free(&lower); return M_NONFINITE; }
            lower.values[row * n + col] = value;
        }
    }
    m_free(&lower);
    return product_finish(product, out);
}

static bool safe_formula_product(double left, double right, double *out) {
    double product = left * right;
    if (!isfinite(product) || (left != 0.0 && right != 0.0 && fabs(product) < DBL_MIN))
        return false;
    *out = product;
    return true;
}

/* Closed forms are useful for small, well-scaled matrices. Fall back when a
 * product loses range or cancellation loses roughly half the significant bits. */
static bool determinant_small(const matrix *source, double *out) {
    const double *v = source->values;
    double terms[6], largest = 0.0, result = 0.0;
    size_t count = source->rows == 2 ? 2 : 6;
    if (source->rows == 2) {
        if (!safe_formula_product(v[0], v[3], &terms[0]) ||
            !safe_formula_product(-v[1], v[2], &terms[1])) return false;
    } else {
        const size_t indices[6][3] = {{0, 4, 8}, {1, 5, 6}, {2, 3, 7},
                                    {2, 4, 6}, {1, 3, 8}, {0, 5, 7}};
        for (size_t i = 0; i < count; i++) {
            double pair;
            if (!safe_formula_product(v[indices[i][0]], v[indices[i][1]], &pair) ||
                !safe_formula_product(pair, v[indices[i][2]], &terms[i])) return false;
            if (i >= 3) terms[i] = -terms[i];
        }
    }
    for (size_t i = 0; i < count; i++) {
        if (fabs(terms[i]) > largest) largest = fabs(terms[i]);
        result += terms[i];
        if (!isfinite(result)) return false;
    }
    if (largest != 0.0 && fabs(result) / largest <= count * sqrt(DBL_EPSILON)) return false;
    if (result != 0.0 && fabs(result) < DBL_MIN) return false;
    *out = result;
    return true;
}

matrix_status m_determinant_with_algorithm(const matrix *source,
                                         matrix_determinant_algorithm algorithm,
                                         double *out) {
    if (!out) return M_ARGUMENT;
    if (algorithm != M_DETERMINANT_AUTO && algorithm != M_DETERMINANT_LU &&
        algorithm != M_DETERMINANT_CHOLESKY) return M_ALGORITHM;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    if (source->rows != source->cols) return M_SHAPE;
    if (algorithm == M_DETERMINANT_LU) return determinant_lu(source, out);
    if (algorithm == M_DETERMINANT_CHOLESKY) return determinant_cholesky(source, out);
    size_t n = source->rows;
    if (n == 0) { *out = 1.0; return M_OK; }
    if (n == 1) { *out = source->values[0]; return M_OK; }
    bool upper = true, lower = true;
    for (size_t row = 0; row < n && (upper || lower); row++)
    for (size_t col = 0; col < row; col++) {
        if (source->values[row * n + col] != 0.0) upper = false;
        if (source->values[col * n + row] != 0.0) lower = false;
    }
    if (upper || lower) {
        determinant_product product = {1.0, 0};
        for (size_t i = 0; i < n; i++) product_append(&product, source->values[i * n + i]);
        return product_finish(product, out);
    }
    if (n <= 3 && determinant_small(source, out)) return M_OK;
    return determinant_lu(source, out);
}

matrix_status m_determinant(const matrix *source, double *out) {
    return m_determinant_with_algorithm(source, M_DETERMINANT_AUTO, out);
}

matrix_status m_triangular(const matrix *source, bool *upper, bool *lower) {
    if (!upper || !lower) return M_ARGUMENT;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    bool u = source->rows == source->cols, l = u;
    if (u) for (size_t row = 0; row < source->rows; row++)
    for (size_t col = 0; col < source->cols; col++) {
        if (row > col && source->values[row * source->cols + col] != 0.0) u = false;
        if (row < col && source->values[row * source->cols + col] != 0.0) l = false;
    }
    *upper = u; *lower = l;
    return M_OK;
}

static double eigen_off_norm(const matrix *work) {
    double norm = 0.0;
    const size_t n = work->rows;
    for (size_t row = 0; row < n; row++)
    for (size_t col = row + 1; col < n; col++) {
        norm = hypot(norm, work->values[row * n + col]);
        norm = hypot(norm, work->values[row * n + col]);
    }
    return norm;
}

matrix_status m_eigen_symmetric(const matrix *source, matrix *values, matrix *vectors) {
    return m_eigen_symmetric_with_options(source, 1e-12, 50, values, vectors);
}

matrix_status m_eigen_symmetric_with_options(const matrix *source, double tolerance,
                                            size_t max_sweeps, matrix *values, matrix *vectors) {
    if (!empty_output(values) || !empty_output(vectors) || values == vectors ||
        source == values || source == vectors || !isfinite(tolerance) ||
        tolerance <= 0.0 || tolerance >= 1.0 || !max_sweeps || max_sweeps > 10000)
        return M_ARGUMENT;
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    if (source->rows != source->cols) return M_SHAPE;
    const size_t n = source->rows;
    double maximum = 0.0;
    bool diagonal = true;
    for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++) {
        const double value = source->values[row * n + col];
        if (value != source->values[col * n + row]) return M_NOT_SYMMETRIC;
        if (row != col && value != 0.0) diagonal = false;
        maximum = fmax(maximum, fabs(value));
    }
    matrix eigenvalues = {0}, eigenvectors = {0}, work = {0};
    status = m_create(n, 1, &eigenvalues);
    if (status != M_OK) goto cleanup;
    status = m_identity(n, &eigenvectors);
    if (status != M_OK) goto cleanup;
    if (diagonal) {
        for (size_t i = 0; i < n; i++) eigenvalues.values[i] = source->values[i * n + i];
    } else {
        status = m_copy(source, &work);
        if (status != M_OK) goto cleanup;
        int exponent = 0;
        (void)frexp(maximum, &exponent);
        double norm = 0.0;
        for (size_t i = 0; i < n * n; i++) {
            work.values[i] = scalbn(work.values[i], -exponent);
            norm = hypot(norm, work.values[i]);
        }
        const double target = tolerance * norm;
        const double skip = target / (2.0 * (double)n);
        size_t sweep = 0;
        while (eigen_off_norm(&work) > target && sweep < max_sweeps) {
            for (size_t p = 0; p < n; p++)
            for (size_t q = p + 1; q < n; q++) {
                const double apq = work.values[p * n + q];
                if (fabs(apq) <= skip) continue;
                const double app = work.values[p * n + p], aqq = work.values[q * n + q];
                const double delta = aqq - app, twice = 2.0 * apq;
                const double t = twice / (delta + copysign(hypot(delta, twice), delta));
                const double c = 1.0 / hypot(1.0, t), s = t * c;
                work.values[p * n + p] = app - t * apq;
                work.values[q * n + q] = aqq + t * apq;
                work.values[p * n + q] = work.values[q * n + p] = 0.0;
                for (size_t k = 0; k < n; k++) {
                    if (k != p && k != q) {
                        const double akp = work.values[k * n + p], akq = work.values[k * n + q];
                        work.values[k * n + p] = work.values[p * n + k] = c * akp - s * akq;
                        work.values[k * n + q] = work.values[q * n + k] = s * akp + c * akq;
                    }
                    const double vkp = eigenvectors.values[k * n + p], vkq = eigenvectors.values[k * n + q];
                    eigenvectors.values[k * n + p] = c * vkp - s * vkq;
                    eigenvectors.values[k * n + q] = s * vkp + c * vkq;
                }
            }
            sweep++;
        }
        if (eigen_off_norm(&work) > target) { status = M_NO_CONVERGENCE; goto cleanup; }
        for (size_t i = 0; i < n; i++) {
            eigenvalues.values[i] = scalbn(work.values[i * n + i], exponent);
            if (!isfinite(eigenvalues.values[i])) { status = M_NONFINITE; goto cleanup; }
        }
    }
    /* Sorting moves whole eigenvectors; signs are otherwise mathematically arbitrary. */
    for (size_t i = 0; i < n; i++) {
        size_t smallest = i;
        for (size_t j = i + 1; j < n; j++)
            if (eigenvalues.values[j] < eigenvalues.values[smallest]) smallest = j;
        if (smallest != i) {
            double swap = eigenvalues.values[i];
            eigenvalues.values[i] = eigenvalues.values[smallest]; eigenvalues.values[smallest] = swap;
            for (size_t row = 0; row < n; row++) {
                swap = eigenvectors.values[row * n + i];
                eigenvectors.values[row * n + i] = eigenvectors.values[row * n + smallest];
                eigenvectors.values[row * n + smallest] = swap;
            }
        }
        size_t largest = 0;
        for (size_t row = 1; row < n; row++)
            if (fabs(eigenvectors.values[row * n + i]) > fabs(eigenvectors.values[largest * n + i])) largest = row;
        if (eigenvectors.values[largest * n + i] < 0.0)
            for (size_t row = 0; row < n; row++) eigenvectors.values[row * n + i] = -eigenvectors.values[row * n + i];
    }
    m_free(&work);
    *values = eigenvalues; *vectors = eigenvectors;
    return M_OK;
cleanup:
    m_free(&work); m_free(&eigenvalues); m_free(&eigenvectors);
    return status;
}

static bool vector3(const matrix *value) {
    return (value->rows == 3 && value->cols == 1) || (value->rows == 1 && value->cols == 3);
}

matrix_status m_cross(const matrix *left, const matrix *right, matrix *out) {
    if (!empty_output(out) || left == out || right == out) return M_ARGUMENT;
    matrix_status status = validate(left);
    if (status != M_OK) return status;
    status = validate(right);
    if (status != M_OK) return status;
    if (!vector3(left) || !vector3(right)) return M_SHAPE;
    matrix result = {0};
    status = m_create(left->rows, left->cols, &result);
    if (status != M_OK) return status;
#ifdef MATRIX_USE_ASM
    m_asm_cross(left->values, right->values, result.values);
#else
    const double *a = left->values, *b = right->values;
    result.values[0] = a[1] * b[2] - a[2] * b[1];
    result.values[1] = a[2] * b[0] - a[0] * b[2];
    result.values[2] = a[0] * b[1] - a[1] * b[0];
#endif
    return finish(&result, out);
}

matrix_status m_rotation_2d(double radians, matrix *out) {
    if (!empty_output(out)) return M_ARGUMENT;
    if (!isfinite(radians)) return M_NONFINITE;
    const double c = cos(radians), s = sin(radians);
    const double entries[] = {c, -s, s, c};
    return m_from_array(2, 2, entries, out);
}

matrix_status m_rotation_x(double radians, matrix *out) {
    if (!empty_output(out)) return M_ARGUMENT;
    if (!isfinite(radians)) return M_NONFINITE;
    const double c = cos(radians), s = sin(radians);
    const double entries[] = {1, 0, 0, 0, c, -s, 0, s, c};
    return m_from_array(3, 3, entries, out);
}

matrix_status m_rotation_y(double radians, matrix *out) {
    if (!empty_output(out)) return M_ARGUMENT;
    if (!isfinite(radians)) return M_NONFINITE;
    const double c = cos(radians), s = sin(radians);
    const double entries[] = {c, 0, s, 0, 1, 0, -s, 0, c};
    return m_from_array(3, 3, entries, out);
}

matrix_status m_rotation_z(double radians, matrix *out) {
    if (!empty_output(out)) return M_ARGUMENT;
    if (!isfinite(radians)) return M_NONFINITE;
    const double c = cos(radians), s = sin(radians);
    const double entries[] = {c, -s, 0, s, c, 0, 0, 0, 1};
    return m_from_array(3, 3, entries, out);
}

matrix_status m_rotation_axis_angle(const matrix *axis, double radians, matrix *out) {
    if (!empty_output(out) || axis == out) return M_ARGUMENT;
    if (!isfinite(radians)) return M_NONFINITE;
    matrix_status status = validate(axis);
    if (status != M_OK) return status;
    if (!vector3(axis)) return M_SHAPE;
    const double maximum = fmax(fabs(axis->values[0]), fmax(fabs(axis->values[1]), fabs(axis->values[2])));
    if (maximum == 0.0) return M_ARGUMENT;
    double x = axis->values[0] / maximum, y = axis->values[1] / maximum, z = axis->values[2] / maximum;
    const double norm = sqrt(x * x + y * y + z * z);
    x /= norm; y /= norm; z /= norm;
    const double c = cos(radians), s = sin(radians), half_sine = sin(radians / 2.0);
    const double t = fabs(radians) < 1.0 ? 2.0 * half_sine * half_sine : 1.0 - c;
    const double entries[] = {
        c + x*x*t, x*y*t - z*s, x*z*t + y*s,
        y*x*t + z*s, c + y*y*t, y*z*t - x*s,
        z*x*t - y*s, z*y*t + x*s, c + z*z*t,
    };
    return m_from_array(3, 3, entries, out);
}


matrix_status m_eigen_general(const matrix *source, matrix *values_real, matrix *values_imag,
                              matrix *vectors_real, matrix *vectors_imag) {
    return m_eigen_general_with_options(source, 1000, values_real, values_imag, vectors_real, vectors_imag);
}

matrix_status m_eigen_general_with_options(const matrix *source, size_t max_iterations,
                                          matrix *values_real, matrix *values_imag,
                                          matrix *vectors_real, matrix *vectors_imag) {
    matrix *outputs[] = {values_real, values_imag, vectors_real, vectors_imag};
    if (!max_iterations || max_iterations > 100000) return M_ARGUMENT;
    for (size_t i = 0; i < 4; i++) {
        if (!empty_output(outputs[i]) || outputs[i] == source) return M_ARGUMENT;
        for (size_t j = 0; j < i; j++) if (outputs[i] == outputs[j]) return M_ARGUMENT;
    }
    matrix_status status = validate(source);
    if (status != M_OK) return status;
    if (source->rows != source->cols) return M_SHAPE;
    const size_t n = source->rows;
    if (n > INT_MAX) return M_MEMORY;
    matrix dr = {0}, di = {0}, vr = {0}, vi = {0}, h = {0}, v = {0}, ort = {0}, powers = {0};
    matrix *buffers[] = {&dr, &di, &vr, &vi, &h, &v, &ort, &powers};
    for (size_t i = 0; i < 8; i++) {
        status = m_create(n, i < 2 || i >= 6 ? 1 : n, buffers[i]);
        if (status != M_OK) goto cleanup;
    }
    switch (la_general_eigen((int)n, source->values, max_iterations, dr.values, di.values,
                            vr.values, vi.values, h.values, v.values, ort.values, powers.values)) {
        case 1: status = M_NO_CONVERGENCE; goto cleanup;
        case 2: status = M_NONFINITE; goto cleanup;
        default: break;
    }
    m_free(&h); m_free(&v); m_free(&ort); m_free(&powers);
    *values_real = dr; *values_imag = di; *vectors_real = vr; *vectors_imag = vi;
    return M_OK;
cleanup:
    for (size_t i = 0; i < 8; i++) m_free(buffers[i]);
    return status;
}

#include "solve_core.h"
struct matrix_factor { la_factor *core; };
static matrix_status solver_status(int status) {
    switch(status){
        case LA_OK:return M_OK;case LA_SHAPE:return M_SHAPE;case LA_MEMORY:return M_MEMORY;
        case LA_RANGE:return M_SOLVER_RANGE;case LA_SINGULAR:return M_SINGULAR;
        case LA_NOT_SPD:return M_NOT_POSITIVE_DEFINITE;case LA_RANK:return M_RANK_DEFICIENT;
        default:return M_ARGUMENT;
    }
}
matrix_status m_factorize(const matrix *source,matrix_factor_algorithm algorithm,matrix_factor **out){
    if(!out||*out)return M_ARGUMENT;
    matrix_status status=validate(source);if(status!=M_OK)return status;
    matrix_factor *factor=calloc(1,sizeof(matrix_factor));if(!factor)return M_MEMORY;
    status=solver_status(la_factor_create(source->rows,source->cols,source->values,algorithm,&factor->core));
    if(status!=M_OK){free(factor);return status;}*out=factor;return M_OK;
}
void m_factor_free(matrix_factor **factor){
    if(factor&&*factor){la_factor_destroy((*factor)->core);free(*factor);*factor=NULL;}
}
matrix_status m_factor_solve(const matrix_factor *factor,const matrix *rhs,matrix *out){
    if(!factor||!factor->core||!empty_output(out)||rhs==out)return M_ARGUMENT;
    matrix_status status=validate(rhs);if(status!=M_OK)return status;
    if(rhs->rows!=factor->core->rows)return M_SHAPE;
    matrix result={0};status=m_create(factor->core->cols,rhs->cols,&result);if(status!=M_OK)return status;
    status=solver_status(la_factor_solve(factor->core,rhs->rows,rhs->cols,rhs->values,result.values,1));
    if(status!=M_OK){m_free(&result);return status;}*out=result;return M_OK;
}
matrix_status m_factor_rcond(const matrix_factor *factor,double *out){
    if(!factor||!out)return M_ARGUMENT;return solver_status(la_factor_rcond(factor->core,out));
}
static matrix_status solve_with(const matrix *source,const matrix *rhs,matrix *out,matrix_factor_algorithm algorithm){
    if(source==out||rhs==out||!empty_output(out))return M_ARGUMENT;
    matrix_factor *factor=NULL;matrix_status status=m_factorize(source,algorithm,&factor);
    if(status==M_OK)status=m_factor_solve(factor,rhs,out);m_factor_free(&factor);return status;
}
matrix_status m_solve(const matrix *source,const matrix *rhs,matrix *out){return solve_with(source,rhs,out,M_LU);}
matrix_status m_least_squares(const matrix *source,const matrix *rhs,matrix *out){return solve_with(source,rhs,out,M_QR);}

#include "svd_core.h"
matrix_status m_svd(const matrix *source,double tolerance,size_t max_sweeps,matrix *u,matrix *values,matrix *vt){
    if(!source||!empty_output(u)||!empty_output(values)||!empty_output(vt)||u==values||u==vt||values==vt||source==u||source==values||source==vt)return M_ARGUMENT;
    matrix_status checked=validate(source);if(checked!=M_OK)return checked;
    size_t k=source->rows<source->cols?source->rows:source->cols;
    matrix a={0},s={0},b={0};
    matrix_status status=m_create(source->rows,k,&a);
    if(status==M_OK)status=m_create(k,1,&s);
    if(status==M_OK)status=m_create(k,source->cols,&b);
    if(status==M_OK){
        int code=la_svd(source->rows,source->cols,source->values,tolerance,max_sweeps,a.values,s.values,b.values);
        status=code==1?M_ARGUMENT:code==2?M_MEMORY:code==3?M_NONFINITE:code==4?M_NO_CONVERGENCE:M_OK;
    }
    if(status==M_OK){*u=a;*values=s;*vt=b;}else{m_free(&a);m_free(&s);m_free(&b);}
    return status;
}

#include "pseudoinverse_core.h"
static double inverse_cutoff(const matrix *a,double cutoff){
    return cutoff==-1?(double)(a->rows>a->cols?a->rows:a->cols)*DBL_EPSILON:cutoff;
}
static matrix_status inverse_apply(const matrix *a,const matrix *rhs,double cutoff,double lambda,matrix *out){
    if(!isfinite(lambda)||lambda<0)return M_ARGUMENT;
    if(!a||!empty_output(out)||a==out||rhs==out)return M_ARGUMENT;
    matrix_status status=validate(a);if(status!=M_OK)return status;
    cutoff=inverse_cutoff(a,cutoff);
    if(!isfinite(cutoff)||cutoff<0||cutoff>1)return M_ARGUMENT;
    if(rhs){status=validate(rhs);if(status!=M_OK)return status;if(rhs->rows!=a->rows)return M_SHAPE;}
    matrix u={0},s={0},vt={0},result={0};
    status=m_svd(a,1e-12,100,&u,&s,&vt);
    size_t cols=rhs?rhs->cols:a->rows;
    if(status==M_OK)status=m_create(a->cols,cols,&result);
    if(status==M_OK){
        int code=la_apply_inverse(a->rows,a->cols,s.rows,u.values,s.values,vt.values,cutoff,rhs?rhs->values:NULL,cols,rhs!=NULL,result.values,lambda);
        if(code)status=code==1?M_ARGUMENT:M_SOLVER_RANGE;
    }
    m_free(&u);m_free(&s);m_free(&vt);
    if(status==M_OK)*out=result;else m_free(&result);
    return status;
}
matrix_status m_pseudoinverse(const matrix *a,double cutoff,matrix *out){return inverse_apply(a,NULL,cutoff,0,out);}
matrix_status m_solve_minimum_norm(const matrix *a,const matrix *b,double cutoff,matrix *out){
    if(!b)return M_ARGUMENT;
    return inverse_apply(a,b,cutoff,0,out);
}
matrix_status m_spectral_diagnostics(const matrix *a,double cutoff,matrix_spectral_diagnostics *out){
    if(!a||!out)return M_ARGUMENT;
    cutoff=inverse_cutoff(a,cutoff);
    if(!isfinite(cutoff)||cutoff<0||cutoff>1)return M_ARGUMENT;
    matrix u={0},s={0},vt={0};matrix_status status=m_svd(a,1e-12,100,&u,&s,&vt);
    if(status==M_OK){
        size_t rank=la_spectral_rank(s.rows,s.values,cutoff);
        matrix_spectral_diagnostics result={rank,s.rows&&s.values[0]?s.values[s.rows-1]/s.values[0]:0,rank?s.values[rank-1]/s.values[0]:0};
        *out=result;
    }
    m_free(&u);m_free(&s);m_free(&vt);return status;
}

matrix_status m_solve_ridge(const matrix *a,const matrix *b,double lambda,matrix *out){
    if(!b)return M_ARGUMENT;
    return inverse_apply(a,b,lambda==0?-1:0,lambda,out);
}

#include "sparse_impl.h"
