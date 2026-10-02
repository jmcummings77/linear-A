#include "trace.h"

_Static_assert(sizeof(m_trace_step) == 8 * sizeof(double), "trace fields must be contiguous");

static m_trace_step steps[M_TRACE_CAPACITY];
static size_t count;
static bool active;
static bool overflow;

void m_trace_clear(void) {
    count = 0;
    active = false;
    overflow = false;
}

void m_trace_begin(void) {
    m_trace_clear();
    active = true;
}

bool m_trace_end(bool success) {
    active = false;
    if (!success || overflow) count = 0;
    return !overflow;
}

void m_trace_record(size_t row, size_t col, size_t k, double left, double right,
                    double product, double sum, size_t source_line) {
    if (!active) return;
    if (count == M_TRACE_CAPACITY) { overflow = true; return; }
    steps[count++] = (m_trace_step){
        (double)row, (double)col, (double)k, left, right, product, sum, (double)source_line
    };
}

size_t m_trace_count(void) { return count; }
const m_trace_step *m_trace_data(void) { return steps; }
