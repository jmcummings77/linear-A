.text
.p2align 2
.globl _compiled_dot
_compiled_dot:
    fmov d0, xzr
    cbz x2, 2f
1:
    ldr d1, [x0], #8
    ldr d2, [x1], #8
    fmul d1, d1, d2
    fadd d0, d0, d1
    subs x2, x2, #1
    b.ne 1b
2:
    ret
