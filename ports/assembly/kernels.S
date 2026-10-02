/* Scalar ARM64 kernels. C owns allocation, validation, and determinant pivoting.
 * Separate multiply and add/sub instructions deliberately avoid fused rounding. */
#if defined(__APPLE__)
#define SYMBOL(name) _##name
#else
#define SYMBOL(name) name
#endif

.text
.p2align 2
.globl SYMBOL(m_asm_add)
SYMBOL(m_asm_add):
    cbz x3, 2f
1:  ldr d0, [x0], #8
    ldr d1, [x1], #8
    fadd d0, d0, d1
    str d0, [x2], #8
    subs x3, x3, #1
    b.ne 1b
2:  ret

.p2align 2
.globl SYMBOL(m_asm_subtract)
SYMBOL(m_asm_subtract):
    cbz x3, 2f
1:  ldr d0, [x0], #8
    ldr d1, [x1], #8
    fsub d0, d0, d1
    str d0, [x2], #8
    subs x3, x3, #1
    b.ne 1b
2:  ret

.p2align 2
.globl SYMBOL(m_asm_scale)
SYMBOL(m_asm_scale):
    cbz x2, 2f
1:  ldr d1, [x0], #8
    fmul d1, d1, d0
    str d1, [x1], #8
    subs x2, x2, #1
    b.ne 1b
2:  ret

.p2align 2
.globl SYMBOL(m_asm_transpose)
SYMBOL(m_asm_transpose):
    cbz x2, 3f
    cbz x3, 3f
    mov x4, #0
1:  mov x5, #0
2:  ldr d0, [x0], #8
    madd x6, x5, x2, x4
    str d0, [x1, x6, lsl #3]
    add x5, x5, #1
    cmp x5, x3
    b.lo 2b
    add x4, x4, #1
    cmp x4, x2
    b.lo 1b
3:  ret

.p2align 2
.globl SYMBOL(m_asm_multiply)
SYMBOL(m_asm_multiply):
    cbz x3, 5f
    cbz x5, 5f
    mov x6, #0
1:  mov x7, #0
2:  fmov d0, xzr
    mov x8, #0
    cbz x4, 4f
3:  madd x9, x6, x4, x8
    ldr d1, [x0, x9, lsl #3]
    madd x10, x8, x5, x7
    ldr d2, [x1, x10, lsl #3]
    fmul d1, d1, d2
    fadd d0, d0, d1
    add x8, x8, #1
    cmp x8, x4
    b.lo 3b
4:  madd x9, x6, x5, x7
    str d0, [x2, x9, lsl #3]
    add x7, x7, #1
    cmp x7, x5
    b.lo 2b
    add x6, x6, #1
    cmp x6, x3
    b.lo 1b
5:  ret

.p2align 2
.globl SYMBOL(m_asm_trace)
SYMBOL(m_asm_trace):
    fmov d0, xzr
    cbz x1, 2f
    add x2, x1, #1
    lsl x2, x2, #3
1:  ldr d1, [x0]
    fadd d0, d0, d1
    add x0, x0, x2
    subs x1, x1, #1
    b.ne 1b
2:  ret

.p2align 2
.globl SYMBOL(m_asm_eliminate)
SYMBOL(m_asm_eliminate):
    cbz x2, 2f
1:  ldr d1, [x1], #8
    ldr d2, [x0]
    fmul d1, d1, d0
    fsub d2, d2, d1
    str d2, [x0], #8
    subs x2, x2, #1
    b.ne 1b
2:  ret

/* x0/x1 point to three inputs; x2 points to three independent outputs. */
.p2align 2
.globl SYMBOL(m_asm_cross)
SYMBOL(m_asm_cross):
    ldp d0, d1, [x0]
    ldr d2, [x0, #16]
    ldp d3, d4, [x1]
    ldr d5, [x1, #16]
    fmul d6, d1, d5
    fmul d7, d2, d4
    fsub d6, d6, d7
    str d6, [x2]
    fmul d6, d2, d3
    fmul d7, d0, d5
    fsub d6, d6, d7
    str d6, [x2, #8]
    fmul d6, d0, d4
    fmul d7, d1, d3
    fsub d6, d6, d7
    str d6, [x2, #16]
    ret

#if !defined(__APPLE__)
.section .note.GNU-stack,"",%progbits
#endif
