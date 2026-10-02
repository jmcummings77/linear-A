#include <cstddef>
extern "C" double compiled_dot(const double *a, const double *b, std::size_t n) {
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) sum += a[i] * b[i];
    return sum;
}
