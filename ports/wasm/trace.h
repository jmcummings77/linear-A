#ifndef LINEAR_A_WASM_TRACE_H
#define LINEAR_A_WASM_TRACE_H

#include <stdbool.h>
#include <stddef.h>

/* Trace builds are single-threaded; capture is bounded and opt-in per call. */
#define M_TRACE_MAX_DIMENSION 8
#define M_TRACE_CAPACITY (M_TRACE_MAX_DIMENSION * M_TRACE_MAX_DIMENSION * M_TRACE_MAX_DIMENSION)

/* All fields are doubles so the JS boundary has one stable, aligned layout. */
typedef struct {
    double row, col, k, left, right, product, sum, source_line;
} m_trace_step;

void m_trace_clear(void);
void m_trace_begin(void);
bool m_trace_end(bool success);
void m_trace_record(size_t row, size_t col, size_t k, double left, double right,
                    double product, double sum, size_t source_line);
size_t m_trace_count(void);
const m_trace_step *m_trace_data(void);

#endif
