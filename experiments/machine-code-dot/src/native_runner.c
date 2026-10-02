#define _DARWIN_C_SOURCE
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <libkern/OSCacheControl.h>
#ifdef RAW_MACHINE_CODE
#include "raw_dot.h"
#endif
typedef double (*dot_fn)(const double *, const double *, size_t);
#ifndef RAW_MACHINE_CODE
extern double compiled_dot(const double *, const double *, size_t);
#endif
static uint64_t nanoseconds(void) {
    struct timespec value;
    clock_gettime(CLOCK_MONOTONIC, &value);
    return (uint64_t)value.tv_sec * 1000000000ULL + (uint64_t)value.tv_nsec;
}
static dot_fn load_dot(void) {
#ifdef RAW_MACHINE_CODE
    void *memory = mmap(NULL, 4096, PROT_READ | PROT_WRITE | PROT_EXEC,
                        MAP_PRIVATE | MAP_ANON | MAP_JIT, -1, 0);
    if (memory == MAP_FAILED) { perror("mmap MAP_JIT"); exit(2); }
    pthread_jit_write_protect_np(0);
    memcpy(memory, raw_dot_words, sizeof(raw_dot_words));
    sys_icache_invalidate(memory, sizeof(raw_dot_words));
    pthread_jit_write_protect_np(1);
    return (dot_fn)memory;
#else
    return compiled_dot;
#endif
}
static size_t number(const char *arg) {
    char *end;
    errno = 0;
    unsigned long value = strtoul(arg, &end, 10);
    if (errno || *end || *arg == '-' || end == arg) exit(2);
    return value;
}
int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: runner length iterations seed\n"); return 2; }
    size_t n = number(argv[1]), iterations = number(argv[2]), seed = number(argv[3]);
    if (n > 1048576 || !iterations || iterations > 100000000 || seed > 1000000) return 2;
    size_t count = (n ? n : 1) * 32;
    double *a = calloc(count, sizeof(double)), *b = calloc(count, sizeof(double));
    if (!a || !b) return 2;
    for (size_t i = 0; i < count; ++i) {
        a[i] = ((int)((i * 17 + seed * 13) % 101) - 50) / 16.0;
        b[i] = ((int)((i * 29 + (seed + 1) * 7) % 103) - 51) / 16.0;
    }
    /* A volatile load keeps all native controls on the same indirect call path. */
    dot_fn volatile selected_dot = load_dot();
    dot_fn dot = selected_dot;
    volatile double sink = 0;
    size_t warm = iterations < 128 ? iterations : 128;
    if (warm < 8) warm = 8;
    for (size_t i = 0; i < warm; ++i) sink = dot(a + (i % 32) * n, b + (i % 32) * n, n);
    double checksum = 0;
    uint64_t start = nanoseconds();
    for (size_t i = 0; i < iterations; ++i) checksum += dot(a + (i % 32) * n, b + (i % 32) * n, n);
    double elapsed = (double)(nanoseconds() - start);
    if (!isfinite(checksum)) return 2;
    printf("{\"elapsed_ns\":%.0f,\"iterations\":%zu,\"checksum\":%.17g}\n", elapsed, iterations, checksum);
    (void)sink;
#ifdef RAW_MACHINE_CODE
    munmap((void *)dot, 4096);
#endif
    free(a); free(b);
    return 0;
}
