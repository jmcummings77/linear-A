// Exercise exactly the same public API and cases as C, with assembly dispatch.
#ifndef MATRIX_USE_ASM
#error Compile this test with -DMATRIX_USE_ASM and kernels.S.
#endif
#include "../c/tests.c"
