#ifndef LINEAR_A_CPP_BRIDGE_HPP
#define LINEAR_A_CPP_BRIDGE_HPP

// Adapter for the shared C/C++ protocol parser and timer. All matrix arithmetic
// delegates to the std::vector-backed C++ API (general QR arithmetic is shared).
#include "matrix.hpp"
#include <new>

enum matrix_status { M_OK, M_ARGUMENT, M_SHAPE, M_INDEX, M_MEMORY, M_NONFINITE, M_ALGORITHM, M_NOT_POSITIVE_DEFINITE, M_NOT_SYMMETRIC, M_NO_CONVERGENCE };
enum matrix_determinant_algorithm { M_DETERMINANT_AUTO, M_DETERMINANT_LU, M_DETERMINANT_CHOLESKY };
struct matrix {
    linear_a::Matrix owned;
    std::size_t rows = 0, cols = 0;
    double *values = nullptr;
    matrix(int = 0) { }
    explicit matrix(linear_a::Matrix value) : owned(std::move(value)) { sync(); }
    matrix(const matrix&) = delete;
    matrix& operator=(const matrix&) = delete;
    matrix& operator=(matrix&& other) noexcept { owned = std::move(other.owned); sync(); return *this; }
    void sync() { rows = owned.rows(); cols = owned.cols(); values = owned.data(); }
};
inline const char *m_error(matrix_status status) {
    switch (status) {
        case M_OK: return "success";
        case M_ARGUMENT: return "invalid matrix argument or incompatible shape";
        case M_SHAPE: return "incompatible matrix dimensions";
        case M_INDEX: return "matrix index out of bounds";
        case M_MEMORY: return "allocation failed or dimensions overflow";
        case M_NONFINITE: return "nonfinite input or arithmetic result";
        case M_ALGORITHM: return "unknown determinant algorithm";
        case M_NOT_POSITIVE_DEFINITE: return "Cholesky requires an exactly symmetric positive definite matrix";
        case M_NOT_SYMMETRIC: return "eigendecomposition requires an exactly symmetric matrix";
        case M_NO_CONVERGENCE: return "eigendecomposition did not converge";
    }
    return "unknown matrix error";
}
template<class Function> matrix_status guarded(Function operation) {
    try { operation(); return M_OK; }
    catch (const std::bad_alloc&) { return M_MEMORY; }
    catch (const std::length_error&) { return M_MEMORY; }
    catch (const std::overflow_error&) { return M_NONFINITE; }
    catch (const std::out_of_range&) { return M_INDEX; }
    catch (const linear_a::NotSymmetricError&) { return M_NOT_SYMMETRIC; }
    catch (const linear_a::EigenConvergenceError&) { return M_NO_CONVERGENCE; }
    catch (const std::domain_error&) { return M_NOT_POSITIVE_DEFINITE; }
    catch (const std::invalid_argument&) { return M_ARGUMENT; }
}
inline matrix_status m_create(std::size_t rows, std::size_t cols, matrix *out) {
    return guarded([&] { *out = matrix(linear_a::Matrix(rows, cols)); });
}
inline void m_free(matrix *value) { *value = matrix(); }
inline matrix_status m_add(const matrix *a, const matrix *b, matrix *out) {
    return guarded([&] { *out = matrix(a->owned + b->owned); });
}
inline matrix_status m_subtract(const matrix *a, const matrix *b, matrix *out) {
    return guarded([&] { *out = matrix(a->owned - b->owned); });
}
inline matrix_status m_scale(const matrix *a, double scalar, matrix *out) {
    return guarded([&] { *out = matrix(a->owned.scale(scalar)); });
}
inline matrix_status m_transpose(const matrix *a, matrix *out) {
    return guarded([&] { *out = matrix(a->owned.transpose()); });
}
inline matrix_status m_multiply(const matrix *a, const matrix *b, matrix *out) {
    return guarded([&] { *out = matrix(a->owned * b->owned); });
}
inline matrix_status m_cross(const matrix *a, const matrix *b, matrix *out) {
    return guarded([&] { *out = matrix(a->owned.cross(b->owned)); });
}
inline matrix_status m_rotation_2d(double radians, matrix *out) {
    return guarded([&] { *out = matrix(linear_a::Matrix::rotation_2d(radians)); });
}
inline matrix_status m_rotation_axis_angle(const matrix *axis, double radians, matrix *out) {
    return guarded([&] { *out = matrix(linear_a::Matrix::rotation_axis_angle(axis->owned, radians)); });
}
inline matrix_status m_trace(const matrix *a, double *value) {
    return guarded([&] { *value = a->owned.trace(); });
}
inline matrix_status m_determinant(const matrix *a, double *value) {
    return guarded([&] { *value = a->owned.determinant(); });
}
inline matrix_status m_determinant_with_algorithm(const matrix *a, matrix_determinant_algorithm algorithm, double *value) {
    linear_a::DeterminantAlgorithm selected;
    switch (algorithm) {
        case M_DETERMINANT_AUTO: selected = linear_a::DeterminantAlgorithm::Auto; break;
        case M_DETERMINANT_LU: selected = linear_a::DeterminantAlgorithm::Lu; break;
        case M_DETERMINANT_CHOLESKY: selected = linear_a::DeterminantAlgorithm::Cholesky; break;
        default: return M_ALGORITHM;
    }
    return guarded([&] { *value = a->owned.determinant(selected); });
}
inline matrix_status m_triangular(const matrix *a, bool *upper, bool *lower) {
    return guarded([&] { auto result = a->owned.triangular(); *upper = result.first; *lower = result.second; });
}
inline matrix_status m_eigen_symmetric(const matrix *a, matrix *values, matrix *vectors) {
    return guarded([&] {
        auto result = a->owned.eigen_symmetric();
        const auto size = result.values.size();
        linear_a::Matrix eigenvalues(size, 1, std::move(result.values));
        *values = matrix(std::move(eigenvalues)); *vectors = matrix(std::move(result.vectors));
    });
}
inline matrix_status m_eigen_general(const matrix *a, matrix *values_real, matrix *values_imag,
                                    matrix *vectors_real, matrix *vectors_imag) {
    return guarded([&] {
        auto result = a->owned.eigen_general();
        const auto n = result.values_real.size();
        *values_real = matrix(linear_a::Matrix(n, 1, std::move(result.values_real)));
        *values_imag = matrix(linear_a::Matrix(n, 1, std::move(result.values_imag)));
        *vectors_real = matrix(std::move(result.vectors_real));
        *vectors_imag = matrix(std::move(result.vectors_imag));
    });
}
enum matrix_factor_algorithm { M_LU=1, M_CHOLESKY=2, M_QR=3 };
using matrix_factor = linear_a::Factorization;
inline matrix_status m_factorize(const matrix *a, matrix_factor_algorithm algorithm, matrix_factor **out) {
    return guarded([&] { *out = new matrix_factor(a->owned, algorithm); });
}
inline void m_factor_free(matrix_factor **factor) { delete *factor; *factor = nullptr; }
inline matrix_status m_factor_solve(const matrix_factor *factor, const matrix *rhs, matrix *out) {
    return guarded([&] { *out = matrix(factor->solve(rhs->owned)); });
}
inline matrix_status m_factor_rcond(const matrix_factor *factor, double *out) {
    return guarded([&] { *out = factor->reciprocal_condition(); });
}
inline matrix_status m_solve(const matrix *a, const matrix *b, matrix *out) {
    return guarded([&] { *out = matrix(a->owned.solve(b->owned)); });
}
inline matrix_status m_least_squares(const matrix *a, const matrix *b, matrix *out) {
    return guarded([&] { *out = matrix(a->owned.least_squares(b->owned)); });
}
#endif
