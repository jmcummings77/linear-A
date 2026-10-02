#define _POSIX_C_SOURCE 200809L
#include "kernels.h"
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void fail(const char *message) { fprintf(stderr,"%s\n",message); exit(1); }
static size_t number(const char *text, size_t maximum) {
    if (!*text || *text == '-') fail("expected a nonnegative integer");
    char *tail; errno=0; unsigned long long value=strtoull(text,&tail,10);
    if (errno || *tail || value>maximum) fail("integer outside experiment limits");
    return (size_t)value;
}
static kernel select_kernel(const char *name) {
    if (!strcmp(name,"ijk")) return multiply_ijk;
    if (!strcmp(name,"ikj")) return multiply_ikj;
    if (!strcmp(name,"blocked")) return multiply_blocked;
#ifdef __aarch64__
    if (!strcmp(name,"neon")) return multiply_neon;
#endif
#ifdef HAVE_ACCELERATE
    if (!strcmp(name,"blas")) return multiply_blas;
#endif
    fail("unsupported kernel"); return NULL;
}
static uint64_t now(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC,&t)) fail("monotonic clock failed");
    return (uint64_t)t.tv_sec*1000000000u+(uint64_t)t.tv_nsec;
}
static double *allocate(size_t count) {
    double *p=malloc((count?count:1)*sizeof(double));
    if (!p) fail("allocation failed");
    return p;
}
static int64_t numerator(size_t index, size_t seed) { return (int64_t)((index*17+seed*13)%33)-16; }
static double checksum(const double *c,size_t count) { return count?c[0]+c[count/2]+c[count-1]:0; }
/* Every entry has an exact integer oracle for the bounded dyadic timing inputs. */
static void verify(const double *c,size_t m,size_t k,size_t n,size_t seed) {
    for (size_t i=0;i<m;i++) for (size_t j=0;j<n;j++) {
        int64_t sum=0;
        for (size_t p=0;p<k;p++) sum+=numerator(i*k+p,seed)*numerator(p*n+j,seed+1);
        if (c[i*n+j] != (double)sum/256) fail("full result differs from exact integer reference");
    }
}
int main(int argc,char **argv) {
    if (argc<6) fail("check KERNEL M K N | bench KERNEL M K N ITERATIONS SEED reused|allocated");
    kernel multiply=select_kernel(argv[2]);
    size_t m=number(argv[3],1024),k=number(argv[4],1024),n=number(argv[5],1024);
    double *a=allocate(m*k),*b=allocate(k*n),*c=allocate(m*n);
    if (!strcmp(argv[1],"check")) {
        if (argc!=6) fail("wrong check argument count");
        for(size_t i=0;i<m*k;i++) if(scanf("%lf",a+i)!=1 || !isfinite(a[i])) fail("expected finite input");
        for(size_t i=0;i<k*n;i++) if(scanf("%lf",b+i)!=1 || !isfinite(b[i])) fail("expected finite input");
        char extra; if(scanf(" %c",&extra)==1) fail("extra input");
        /* Poison output to catch missing initialization, including zero inner dimension. */
        for(size_t i=0;i<m*n;i++) c[i]=NAN;
        multiply(a,b,c,m,k,n);
        for(size_t i=0;i<m*n;i++) if(!isfinite(c[i])) fail("nonfinite output");
        printf("{\"rows\":%zu,\"cols\":%zu,\"values\":[",m,n);
        for(size_t i=0;i<m*n;i++) printf("%s%.17g",i?",":"",c[i]);
        puts("]}");
    } else if (!strcmp(argv[1],"bench")) {
        if(argc!=9) fail("wrong bench argument count");
        size_t iterations=number(argv[6],10000000),seed=number(argv[7],1000000);
        int allocated=!strcmp(argv[8],"allocated");
        if(!iterations || (!allocated && strcmp(argv[8],"reused"))) fail("invalid timing mode or iterations");
        for(size_t i=0;i<m*k;i++) a[i]=numerator(i,seed)/16.0;
        for(size_t i=0;i<k*n;i++) b[i]=numerator(i,seed+1)/16.0;
        multiply(a,b,c,m,k,n); verify(c,m,k,n,seed);
        double expected=checksum(c,m*n); volatile double consumed=0;
        for(size_t i=0;i<3;i++) multiply(a,b,c,m,k,n);
        uint64_t start=now();
        for(size_t i=0;i<iterations;i++) {
            double *out=allocated?allocate(m*n):c;
            multiply(a,b,out,m,k,n);
            consumed+=checksum(out,m*n);
            if(allocated) free(out);
        }
        uint64_t elapsed=now()-start;
        /* Repeat full verification outside timing, using the same buffers. */
        multiply(a,b,c,m,k,n); verify(c,m,k,n,seed);
        if(consumed!=expected*(double)iterations) fail("timed checksum mismatch");
        printf("{\"elapsed_ns\":%llu,\"iterations\":%zu,\"checksum\":%.17g,\"full_result_verified\":true,\"output_allocations\":%zu,\"requested_output_bytes\":%zu}\n",
            (unsigned long long)elapsed,iterations,(double)consumed,allocated?iterations:0,allocated?iterations*m*n*sizeof(double):0);
    } else fail("unknown mode");
    free(a);free(b);free(c);return 0;
}
