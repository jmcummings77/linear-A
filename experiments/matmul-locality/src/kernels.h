#ifndef MATMUL_LOCALITY_H
#define MATMUL_LOCALITY_H
#include <stddef.h>
typedef void (*kernel)(const double *, const double *, double *, size_t, size_t, size_t);
void multiply_ijk(const double *, const double *, double *, size_t, size_t, size_t);
void multiply_ikj(const double *, const double *, double *, size_t, size_t, size_t);
void multiply_blocked(const double *, const double *, double *, size_t, size_t, size_t);
#ifdef __aarch64__
void multiply_neon(const double *, const double *, double *, size_t, size_t, size_t);
#endif
#ifdef HAVE_ACCELERATE
void multiply_blas(const double *, const double *, double *, size_t, size_t, size_t);
#endif
#endif
