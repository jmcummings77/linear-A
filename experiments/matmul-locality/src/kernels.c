/* Unchecked kernels: distinct, correctly sized row-major buffers; finite inputs.
 * Each handwritten output reduces k in ascending order, with no fused multiply-add.
 * Output initialization is part of every kernel. No kernel allocates memory. */
#include "kernels.h"
#include <string.h>
#ifdef __aarch64__
#include <arm_neon.h>
#endif
#ifdef HAVE_ACCELERATE
#include <Accelerate/Accelerate.h>
#endif
#define TILE 32
static size_t end(size_t start, size_t limit) { return start + TILE < limit ? start + TILE : limit; }

void multiply_ijk(const double *restrict a, const double *restrict b, double *restrict c, size_t m, size_t k, size_t n) {
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++) {
            double sum = 0;
            for (size_t p = 0; p < k; p++) sum += a[i*k+p] * b[p*n+j];
            c[i*n+j] = sum;
        }
}

void multiply_ikj(const double *restrict a, const double *restrict b, double *restrict c, size_t m, size_t k, size_t n) {
    memset(c, 0, m*n*sizeof(double));
    for (size_t i = 0; i < m; i++)
        for (size_t p = 0; p < k; p++) {
            double value = a[i*k+p];
            for (size_t j = 0; j < n; j++) c[i*n+j] += value * b[p*n+j];
        }
}

void multiply_blocked(const double *restrict a, const double *restrict b, double *restrict c, size_t m, size_t k, size_t n) {
    memset(c, 0, m*n*sizeof(double));
    for (size_t ii = 0; ii < m; ii += TILE)
        for (size_t pp = 0; pp < k; pp += TILE)
            for (size_t jj = 0; jj < n; jj += TILE)
                for (size_t i = ii; i < end(ii,m); i++)
                    for (size_t p = pp; p < end(pp,k); p++) {
                        double value = a[i*k+p];
                        for (size_t j = jj; j < end(jj,n); j++) c[i*n+j] += value * b[p*n+j];
                    }
}

#ifdef __aarch64__
void multiply_neon(const double *restrict a, const double *restrict b, double *restrict c, size_t m, size_t k, size_t n) {
    memset(c, 0, m*n*sizeof(double));
    for (size_t ii = 0; ii < m; ii += TILE)
        for (size_t pp = 0; pp < k; pp += TILE)
            for (size_t jj = 0; jj < n; jj += TILE)
                for (size_t i = ii; i < end(ii,m); i++)
                    for (size_t p = pp; p < end(pp,k); p++) {
                        const double value = a[i*k+p];
                        const float64x2_t av = vdupq_n_f64(value);
                        size_t j = jj, stop = end(jj,n);
                        for (; j + 1 < stop; j += 2) {
                            float64x2_t cv = vld1q_f64(c+i*n+j);
                            float64x2_t bv = vld1q_f64(b+p*n+j);
                            vst1q_f64(c+i*n+j, vaddq_f64(cv, vmulq_f64(av,bv)));
                        }
                        for (; j < stop; j++) c[i*n+j] += value * b[p*n+j];
                    }
}
#endif
#ifdef HAVE_ACCELERATE
void multiply_blas(const double *restrict a, const double *restrict b, double *restrict c, size_t m, size_t k, size_t n) {
    if (!m || !n) return;
    if (!k) { memset(c,0,m*n*sizeof(double)); return; }
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)m,(int)n,(int)k,
                1,a,(int)k,b,(int)n,0,c,(int)n);
}
#endif
