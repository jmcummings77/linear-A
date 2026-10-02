/* Hand-encoded AArch64 instructions. No compiler/assembler generated these words.
   double dot(const double *a [x0], const double *b [x1], size_t n [x2]);
   Result: d0. Caller-saved x0/x1/x2/d0/d1/d2 and condition flags only. */
static const unsigned int raw_dot_words[] = {
    0x9e6703e0, /* 00: fmov d0, xzr             -- sum = +0 */
    0xb40000e2, /* 04: cbz x2, +28             -- empty -> ret */
    0xfc408401, /* 08: ldr d1, [x0], #8        -- load a; advance */
    0xfc408422, /* 0c: ldr d2, [x1], #8        -- load b; advance */
    0x1e620821, /* 10: fmul d1, d1, d2         -- separate multiply */
    0x1e612800, /* 14: fadd d0, d0, d1         -- left-to-right sum */
    0xf1000442, /* 18: subs x2, x2, #1         -- remaining-- */
    0x54ffff61, /* 1c: b.ne -20                -- back to 08 */
    0xd65f03c0  /* 20: ret */
};
