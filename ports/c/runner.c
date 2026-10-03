#ifndef __APPLE__
#define _POSIX_C_SOURCE 200809L
#endif
#ifdef __cplusplus
#include "../cpp/bridge.hpp"
#else
#include "matrix.h"
#endif
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void fail(const char *message) { fprintf(stderr, "%s\n", message); exit(1); }
static void require(matrix_status status) { if (status != M_OK) fail(m_error(status)); }
static size_t integer(const char *text) {
    if (!text[0] || text[0] == '-') fail("expected a nonnegative integer");
    errno = 0;
    char *end;
    uintmax_t value = strtoumax(text, &end, 10);
    if (errno || *end || end == text || value > SIZE_MAX) fail("invalid integer or dimension overflow");
    return (size_t)value;
}
static double number(const char *text) {
    errno = 0;
    char *end;
    double value = strtod(text, &end);
    if (*end || end == text || !isfinite(value)) fail("expected a finite number");
    return value;
}
static uint64_t now_ns(void) {
#ifdef __APPLE__
    return clock_gettime_nsec_np(CLOCK_UPTIME_RAW);
#else
    struct timespec value;
    if (clock_gettime(CLOCK_MONOTONIC, &value)) fail("monotonic timer failed");
    return (uint64_t)value.tv_sec * UINT64_C(1000000000) + (uint64_t)value.tv_nsec;
#endif
}
static bool binary(const char *op) {
    return !strcmp(op, "add") || !strcmp(op, "subtract") || !strcmp(op, "multiply") || !strcmp(op, "cross") || !strcmp(op,"solve") || !strcmp(op,"solve_cholesky") || !strcmp(op,"least_squares");
}
static bool rotation_operation(const char *op) { return !strcmp(op, "rotation2d") || !strcmp(op, "rotation3d"); }
static bool determinant_operation(const char *op) {
    return !strcmp(op, "determinant") || !strcmp(op, "determinant_lu") ||
           !strcmp(op, "determinant_cholesky") || !strcmp(op, "determinant_spd_lu");
}
static bool spd_operation(const char *op) { return !strcmp(op, "determinant_cholesky") || !strcmp(op, "determinant_spd_lu"); }
static bool scalar_result(const char *op) { return !strcmp(op, "trace") || determinant_operation(op) || !strcmp(op,"rcond"); }
static void read_values(matrix *value) {
    for (size_t i = 0; i < value->rows * value->cols; i++) {
        if (scanf("%lf", &value->values[i]) != 1 || !isfinite(value->values[i])) fail("expected exactly the requested finite matrix values");
    }
}
static void execute(const char *op, const matrix *a, const matrix *b, double scalar, matrix *result, double *value) {
    if (!strcmp(op, "add")) require(m_add(a, b, result));
    else if (!strcmp(op, "subtract")) require(m_subtract(a, b, result));
    else if (!strcmp(op, "multiply")) require(m_multiply(a, b, result));
    else if (!strcmp(op,"svd") || !strcmp(op,"svd_one_sweep")) {
        matrix u={0},s={0},vt={0};
        require(m_svd(a,1e-12,!strcmp(op,"svd")?100:1,&u,&s,&vt));
        size_t k=s.rows;
        require(m_create(a->rows+1+a->cols,k,result));
        for(size_t i=0;i<a->rows*k;i++)result->values[i]=u.values[i];
        for(size_t j=0;j<k;j++)result->values[a->rows*k+j]=s.values[j];
        for(size_t i=0;i<a->cols;i++)for(size_t j=0;j<k;j++)result->values[(a->rows+1+i)*k+j]=vt.values[j*a->cols+i];
        m_free(&u);m_free(&s);m_free(&vt);
    }
    else if (!strcmp(op, "solve")) require(m_solve(a,b,result));
    else if (!strcmp(op, "least_squares")) require(m_least_squares(a,b,result));
    else if (!strcmp(op, "solve_cholesky") || !strcmp(op,"rcond")) {
        matrix_factor *factor = NULL;
        require(m_factorize(a, !strcmp(op,"rcond") ? M_LU : M_CHOLESKY, &factor));
        matrix_status status = !strcmp(op,"rcond") ? m_factor_rcond(factor,value) : m_factor_solve(factor,b,result);
        m_factor_free(&factor); require(status);
    }
    else if (!strcmp(op, "cross")) require(m_cross(a, b, result));
    else if (!strcmp(op, "rotation2d")) {
        if (a->rows || a->cols) fail("rotation2d requires empty 0 by 0 input");
        require(m_rotation_2d(scalar, result));
    }
    else if (!strcmp(op, "rotation3d")) require(m_rotation_axis_angle(a, scalar, result));
    else if (!strcmp(op, "scale")) require(m_scale(a, scalar, result));
    else if (!strcmp(op, "transpose")) require(m_transpose(a, result));
    else if (!strcmp(op, "trace")) require(m_trace(a, value));
    else if (!strcmp(op, "determinant")) require(m_determinant(a, value));
    else if (!strcmp(op, "determinant_lu") || !strcmp(op, "determinant_spd_lu")) require(m_determinant_with_algorithm(a, M_DETERMINANT_LU, value));
    else if (!strcmp(op, "determinant_cholesky")) require(m_determinant_with_algorithm(a, M_DETERMINANT_CHOLESKY, value));
    else fail("unknown operation");
}
static void print_matrix(const matrix *value) {
    printf("{\"rows\":%zu,\"cols\":%zu,\"values\":[", value->rows, value->cols);
    for (size_t i = 0; i < value->rows * value->cols; i++) printf("%s%.17g", i ? "," : "", value->values[i]);
    printf("]}");
}
static int check(int argc, char **argv) {
    if (argc < 5) fail("usage: runner check OP ROWS COLS [BROWS BCOLS | SCALAR]");
    const char *op = argv[2];
    int expected = binary(op) ? 7 : (!strcmp(op, "scale") || rotation_operation(op)) ? 6 : 5;
    if (argc != expected) fail("incorrect argument count for operation");
    matrix a = {0}, b = {0}, result = {0};
    require(m_create(integer(argv[3]), integer(argv[4]), &a));
    if (binary(op)) require(m_create(integer(argv[5]), integer(argv[6]), &b));
    double scalar = (!strcmp(op, "scale") || rotation_operation(op)) ? number(argv[5]) : 0.0;
    read_values(&a); read_values(&b);
    int trailing;
    do { trailing = getchar(); } while (trailing != EOF && isspace((unsigned char)trailing));
    if (trailing != EOF) fail("unexpected extra input");
    if (!strcmp(op, "triangular")) {
        bool upper, lower;
        require(m_triangular(&a, &upper, &lower));
        printf("{\"upper\":%s,\"lower\":%s}\n", upper ? "true" : "false", lower ? "true" : "false");
    } else if (!strcmp(op, "eigen_general")) {
        matrix imag = {0}, vr = {0}, vi = {0};
        require(m_eigen_general(&a, &result, &imag, &vr, &vi));
        printf("{\"eigenvalues_real\":[");
        for (size_t i = 0; i < result.rows; i++) printf("%s%.17g", i ? "," : "", result.values[i]);
        printf("],\"eigenvalues_imag\":[");
        for (size_t i = 0; i < imag.rows; i++) printf("%s%.17g", i ? "," : "", imag.values[i]);
        printf("],\"eigenvectors_real\":"); print_matrix(&vr);
        printf(",\"eigenvectors_imag\":"); print_matrix(&vi); puts("}");
        m_free(&imag); m_free(&vr); m_free(&vi);
    } else if (!strcmp(op, "eigen_symmetric")) {
        matrix vectors = {0};
        require(m_eigen_symmetric(&a, &result, &vectors));
        printf("{\"eigenvalues\":[");
        for (size_t i = 0; i < result.rows; i++) printf("%s%.17g", i ? "," : "", result.values[i]);
        printf("],\"eigenvectors\":");
        printf("{\"rows\":%zu,\"cols\":%zu,\"values\":[", vectors.rows, vectors.cols);
        for (size_t i = 0; i < vectors.rows * vectors.cols; i++) printf("%s%.17g", i ? "," : "", vectors.values[i]);
        puts("]}}");
        m_free(&vectors);
    } else {
        double value = 0.0;
        execute(op, &a, &b, scalar, &result, &value);
        if (scalar_result(op)) printf("{\"value\":%.17g}\n", value);
        else { print_matrix(&result); puts(""); }
    }
    m_free(&result); m_free(&a); m_free(&b);
    return 0;
}
static double one_iteration(const char *op, const matrix *a, const matrix *b) {
    // Keep repeated scalar operations from being hoisted out of the timed loop.
    // This is a compiler barrier, with no runtime instruction on GCC/Clang.
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : : "r"(a), "r"(b) : "memory");
#endif
    matrix result = {0};
    double value = 0.0;
    if (!strcmp(op, "eigen_general")) {
        matrix imag = {0}, vr = {0}, vi = {0};
        require(m_eigen_general(a, &result, &imag, &vr, &vi));
        for (size_t i = 0; i < result.rows; i++) value += (double)(i + 1) * (result.values[i] + fabs(imag.values[i]));
        for (size_t i = 0; i < vr.rows * vr.cols; i++) value += vr.values[i] * vr.values[i] + vi.values[i] * vi.values[i];
        m_free(&result); m_free(&imag); m_free(&vr); m_free(&vi);
        if (!isfinite(value)) fail("nonfinite checksum");
        return value;
    }
    if (!strcmp(op, "eigen_symmetric")) {
        matrix vectors = {0};
        require(m_eigen_symmetric(a, &result, &vectors));
        for (size_t i = 0; i < result.rows; i++) value += (double)(i + 1) * result.values[i];
        for (size_t i = 0; i < vectors.rows * vectors.cols; i++) value += vectors.values[i] * vectors.values[i];
        m_free(&result); m_free(&vectors);
        if (!isfinite(value)) fail("nonfinite checksum");
        return value;
    }
    execute(op, a, b, rotation_operation(op) ? 0.5 : 1.25, &result, &value);
    if (!scalar_result(op)) {
        size_t count = result.rows * result.cols;
        value = count ? result.values[0] + result.values[count / 2] + result.values[count - 1] : 0.0;
    }
    m_free(&result);
    if (!isfinite(value)) fail("nonfinite checksum");
    return value;
}
static int bench(int argc, char **argv) {
    if (argc != 6) fail("usage: runner bench OP SIZE ITERATIONS SEED");
    const char *op = argv[2];
    size_t size = integer(argv[3]), iterations = integer(argv[4]), seed = integer(argv[5]);
    if (!size || !iterations || seed > 2147483646U) fail("positive size and iterations, seed 0..2147483646 required");
    if ((!strcmp(op, "cross") || !strcmp(op, "rotation3d")) && size != 3) fail("cross and rotation3d require benchmark size 3");
    if (!strcmp(op, "rotation2d") && size != 2) fail("rotation2d requires benchmark size 2");
    matrix a = {0}, b = {0};
    const size_t rows = !strcmp(op, "rotation2d") ? 0 : size;
    const size_t cols = !strcmp(op, "rotation2d") ? 0 : (!strcmp(op, "cross") || !strcmp(op, "rotation3d")) ? 1 : size;
    require(m_create(rows, cols, &a)); require(m_create(rows, cols, &b));
    for (size_t i = 0; i < rows * cols; i++) {
        a.values[i] = ((double)(((i % 101) * 17 + (seed % 101) * 13) % 101) - 50.0) / 16.0;
        b.values[i] = ((double)(((i % 101) * 17 + ((seed + 1) % 101) * 13) % 101) - 50.0) / 16.0;
    }
    if (!strcmp(op, "rotation3d")) { a.values[0] = 1; a.values[1] = 2; a.values[2] = 3; }
    if (!strcmp(op, "cross")) b.values[2] = -b.values[2];
    if (spd_operation(op)) for (size_t row = 0; row < size; row++)
    for (size_t col = 0; col < size; col++) {
        double forward = (double)((((row * size + col) % 101) * 17 + (seed % 101) * 13) % 101) - 50.0;
        double reverse = (double)((((col * size + row) % 101) * 17 + (seed % 101) * 13) % 101) - 50.0;
        a.values[row * size + col] = (forward + reverse) / 32.0;
    }
    if (determinant_operation(op)) for (size_t i = 0; i < size; i++) a.values[i * size + i] += (double)size * 4.0;
    if (!strcmp(op, "eigen_symmetric")) for (size_t row = 0; row < size; row++)
    for (size_t col = 0; col < size; col++)
        a.values[row * size + col] = row == col ? 2.0 + (double)(seed % 17) / 16.0 :
            (row + 1 == col || col + 1 == row) ? -1.0 : 0.0;
    if (!strcmp(op, "eigen_general")) for (size_t row = 0; row < size; row++)
    for (size_t col = 0; col < size; col++) {
        const double block = (double)(row / 2);
        a.values[row * size + col] = row == col ? 1.0 + (double)(seed % 17) / 16.0 + block / 8.0 :
            (!(row % 2) && col == row + 1) ? -(0.5 + block / 16.0) :
            ((row % 2) && col + 1 == row) ? 0.5 + block / 16.0 :
            row < col ? ((double)((row * 3 + col * 5 + seed) % 11) - 5.0) / 32.0 : 0.0;
    }
    size_t warmup = iterations < 100 ? iterations : 100;
    if (warmup < 5) warmup = 5;
    for (size_t i = 0; i < warmup; i++) (void)one_iteration(op, &a, &b);
    double checksum = 0.0;
    uint64_t start = now_ns();
    for (size_t i = 0; i < iterations; i++) checksum += one_iteration(op, &a, &b);
    uint64_t elapsed = now_ns() - start;
    if (!isfinite(checksum)) fail("nonfinite accumulated checksum");
    printf("{\"elapsed_ns\":%" PRIu64 ",\"iterations\":%zu,\"checksum\":%.17g}\n", elapsed, iterations, checksum);
    m_free(&a); m_free(&b);
    return 0;
}
int main(int argc, char **argv) {
    if (argc < 2) fail("expected check or bench command");
    if (!strcmp(argv[1], "check")) return check(argc, argv);
    if (!strcmp(argv[1], "bench")) return bench(argc, argv);
    fail("expected check or bench command");
    return 1;
}
