#ifndef LINEAR_A_CPP_MATRIX_HPP
#define LINEAR_A_CPP_MATRIX_HPP

#include <algorithm>
#include <memory>
#include "../c/solve_core.h"
#include "../c/svd_core.h"
#include "../c/pseudoinverse_core.h"
#include "../c/general_eigen.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace linear_a {

enum class DeterminantAlgorithm { Auto, Lu, Cholesky };
class Factorization;
struct SymmetricEigenResult;
struct GeneralEigenResult;
class NotSymmetricError : public std::invalid_argument {
public: NotSymmetricError() : std::invalid_argument("eigendecomposition requires an exactly symmetric matrix") { }
};
class EigenConvergenceError : public std::runtime_error {
public: EigenConvergenceError() : std::runtime_error("eigendecomposition did not converge") { }
};

// Owned row-major float64 matrices. Results never mutate their operands.
// Arithmetic rejects nonfinite results. Determinants use partial pivoting and
// retain ordinary float64 rounding/underflow limitations, not exact arithmetic.
class Matrix {
    std::size_t rows_ = 0, cols_ = 0;
    std::vector<double> values_;

    static std::size_t count(std::size_t rows, std::size_t cols) {
        if (cols && rows > std::numeric_limits<std::size_t>::max() / cols / sizeof(double))
            throw std::length_error("matrix dimensions overflow");
        return rows * cols;
    }
    static double finite(double value) {
        if (!std::isfinite(value)) throw std::overflow_error("nonfinite input or arithmetic result");
        return value;
    }
    void validate() const { for (double value : values_) finite(value); }
    void square() const { if (rows_ != cols_) throw std::invalid_argument("matrix must be square"); }
    Matrix combine(const Matrix& other, bool subtract) const {
        validate(); other.validate();
        if (rows_ != other.rows_ || cols_ != other.cols_) throw std::invalid_argument("matrix dimensions must match");
        Matrix result(rows_, cols_);
        for (std::size_t i = 0; i < values_.size(); i++)
            result.values_[i] = subtract ? values_[i] - other.values_[i] : values_[i] + other.values_[i];
        result.validate();
        return result;
    }
    struct DeterminantProduct {
        double fraction = 1.0;
        std::int64_t exponent = 0;
        void append(double factor) {
            int factor_exponent = 0, normalization = 0;
            fraction = std::frexp(fraction * std::frexp(factor, &factor_exponent), &normalization);
            exponent += factor_exponent + normalization;
        }
        double finish() const {
            if (fraction == 0.0) return fraction;
            if (exponent > std::numeric_limits<double>::max_exponent)
                throw std::overflow_error("nonfinite determinant result");
            if (exponent < std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits)
                return std::copysign(0.0, fraction);
            return finite(std::scalbn(fraction, static_cast<int>(exponent)));
        }
    };
    static double quotient_product(double numerator, double denominator, double value) {
        int numerator_exponent, denominator_exponent, value_exponent;
        const double fraction = std::frexp(numerator, &numerator_exponent) /
                                std::frexp(denominator, &denominator_exponent) *
                                std::frexp(value, &value_exponent);
        return std::scalbn(fraction, numerator_exponent - denominator_exponent + value_exponent);
    }
    double determinant_lu() const {
        auto work = values_;
        DeterminantProduct product;
        const std::size_t n = rows_;
        for (std::size_t k = 0; k < n; k++) {
            std::size_t pivot = k;
            for (std::size_t row = k + 1; row < n; row++)
                if (std::abs(work[row * n + k]) > std::abs(work[pivot * n + k])) pivot = row;
            if (work[pivot * n + k] == 0.0) return 0.0;
            if (pivot != k) {
                for (std::size_t col = 0; col < n; col++) std::swap(work[k * n + col], work[pivot * n + col]);
                product.fraction = -product.fraction;
            }
            const double diagonal = work[k * n + k];
            product.append(diagonal);
            for (std::size_t row = k + 1; row < n; row++) {
                const double numerator = work[row * n + k], factor = numerator / diagonal;
                work[row * n + k] = 0.0;
                // Delay division range loss until after the pivot-row product.
                const bool tiny_factor = numerator != 0.0 && std::abs(factor) < std::numeric_limits<double>::min();
                for (std::size_t col = k + 1; col < n; col++) {
                    const double update = tiny_factor ? quotient_product(numerator, diagonal, work[k * n + col])
                                                      : factor * work[k * n + col];
                    work[row * n + col] = finite(work[row * n + col] - update);
                }
            }
        }
        return product.finish();
    }
    double determinant_cholesky() const {
        const std::size_t n = rows_;
        const auto reject = [] { throw std::domain_error("Cholesky requires an exactly symmetric positive definite matrix"); };
        for (std::size_t row = 0; row < n; row++)
        for (std::size_t col = 0; col < row; col++)
            if (values_[row * n + col] != values_[col * n + row]) reject();
        std::vector<double> lower(values_.size(), 0.0);
        DeterminantProduct product;
        for (std::size_t row = 0; row < n; row++)
        for (std::size_t col = 0; col <= row; col++) {
            double value = values_[row * n + col];
            for (std::size_t k = 0; k < col; k++) value -= lower[row * n + k] * lower[col * n + k];
            finite(value);
            if (row == col) {
                if (value <= 0.0) reject();
                lower[row * n + col] = std::sqrt(value);
                product.append(value); // Squared Cholesky diagonal, before sqrt rounding.
            } else lower[row * n + col] = finite(value / lower[col * n + col]);
        }
        return product.finish();
    }
    static bool safe_formula_product(double left, double right, double& out) {
        const double result = left * right;
        if (!std::isfinite(result) || (left != 0.0 && right != 0.0 && std::abs(result) < std::numeric_limits<double>::min()))
            return false;
        out = result;
        return true;
    }
    bool determinant_small(double& out) const {
        double terms[6], result = 0.0, largest = 0.0;
        const std::size_t count = rows_ == 2 ? 2 : 6;
        if (rows_ == 2) {
            if (!safe_formula_product(values_[0], values_[3], terms[0]) ||
                !safe_formula_product(-values_[1], values_[2], terms[1])) return false;
        } else {
            constexpr std::size_t indices[6][3] = {{0, 4, 8}, {1, 5, 6}, {2, 3, 7},
                                                  {2, 4, 6}, {1, 3, 8}, {0, 5, 7}};
            for (std::size_t i = 0; i < count; i++) {
                double pair;
                if (!safe_formula_product(values_[indices[i][0]], values_[indices[i][1]], pair) ||
                    !safe_formula_product(pair, values_[indices[i][2]], terms[i])) return false;
                if (i >= 3) terms[i] = -terms[i];
            }
        }
        for (std::size_t i = 0; i < count; i++) {
            largest = std::max(largest, std::abs(terms[i]));
            result += terms[i];
            if (!std::isfinite(result)) return false;
        }
        if (largest != 0.0 && std::abs(result) / largest <= count * std::sqrt(std::numeric_limits<double>::epsilon())) return false;
        if (result != 0.0 && std::abs(result) < std::numeric_limits<double>::min()) return false;
        out = result;
        return true;
    }
public:
    Matrix() = default;
    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default;
    Matrix(Matrix&& other) noexcept
        : rows_(std::exchange(other.rows_, 0)), cols_(std::exchange(other.cols_, 0)), values_(std::move(other.values_)) { }
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            rows_ = std::exchange(other.rows_, 0);
            cols_ = std::exchange(other.cols_, 0);
            values_ = std::move(other.values_);
        }
        return *this;
    }
    Matrix(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols), values_(count(rows, cols), 0.0) { }
    Matrix(std::size_t rows, std::size_t cols, std::vector<double> values)
        : rows_(rows), cols_(cols), values_(std::move(values)) {
        if (values_.size() != count(rows, cols)) throw std::invalid_argument("value count does not match dimensions");
        validate();
    }
    static Matrix identity(std::size_t size) {
        Matrix result(size, size);
        for (std::size_t i = 0; i < size; i++) result.values_[i * size + i] = 1.0;
        return result;
    }
    // Right-handed active rotations of column vectors, with angles in radians.
    static Matrix rotation_2d(double radians) {
        finite(radians);
        const double c = std::cos(radians), s = std::sin(radians);
        return Matrix(2, 2, {c, -s, s, c});
    }
    static Matrix rotation_x(double radians) {
        finite(radians);
        const double c = std::cos(radians), s = std::sin(radians);
        return Matrix(3, 3, {1,0,0, 0,c,-s, 0,s,c});
    }
    static Matrix rotation_y(double radians) {
        finite(radians);
        const double c = std::cos(radians), s = std::sin(radians);
        return Matrix(3, 3, {c,0,s, 0,1,0, -s,0,c});
    }
    static Matrix rotation_z(double radians) {
        finite(radians);
        const double c = std::cos(radians), s = std::sin(radians);
        return Matrix(3, 3, {c,-s,0, s,c,0, 0,0,1});
    }
    static Matrix rotation_axis_angle(const Matrix& axis, double radians) {
        finite(radians); axis.validate();
        if (!((axis.rows_ == 3 && axis.cols_ == 1) || (axis.rows_ == 1 && axis.cols_ == 3)))
            throw std::invalid_argument("rotation axis must be a three-dimensional vector");
        const double maximum = std::max({std::abs(axis.values_[0]), std::abs(axis.values_[1]), std::abs(axis.values_[2])});
        if (maximum == 0.0) throw std::invalid_argument("rotation axis must be nonzero");
        double x = axis.values_[0] / maximum, y = axis.values_[1] / maximum, z = axis.values_[2] / maximum;
        const double norm = std::sqrt(x*x + y*y + z*z);
        x /= norm; y /= norm; z /= norm;
        const double c = std::cos(radians), s = std::sin(radians), half_sine = std::sin(radians / 2.0);
        const double t = std::abs(radians) < 1.0 ? 2.0 * half_sine * half_sine : 1.0 - c;
        return Matrix(3, 3, {c+x*x*t, x*y*t-z*s, x*z*t+y*s,
                             y*x*t+z*s, c+y*y*t, y*z*t-x*s,
                             z*x*t-y*s, z*y*t+x*s, c+z*z*t});
    }
    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    const std::vector<double>& values() const noexcept { return values_; }
    // Useful for bulk interoperation. Callers preserve the vector's size; every
    // arithmetic operation validates its stored values before calculating.
    double* data() noexcept { return values_.data(); }
    double get(std::size_t row, std::size_t col) const {
        if (row >= rows_ || col >= cols_) throw std::out_of_range("matrix index out of bounds");
        return values_[row * cols_ + col];
    }
    void set(std::size_t row, std::size_t col, double value) {
        if (row >= rows_ || col >= cols_) throw std::out_of_range("matrix index out of bounds");
        values_[row * cols_ + col] = finite(value);
    }
    std::vector<double> row(std::size_t index) const {
        if (index >= rows_) throw std::out_of_range("row index out of bounds");
        auto start = values_.begin() + index * cols_;
        return {start, start + cols_};
    }
    std::vector<double> column(std::size_t index) const {
        if (index >= cols_) throw std::out_of_range("column index out of bounds");
        std::vector<double> result(rows_);
        for (std::size_t row = 0; row < rows_; row++) result[row] = values_[row * cols_ + index];
        return result;
    }
    Matrix operator+(const Matrix& other) const { return combine(other, false); }
    Matrix operator-(const Matrix& other) const { return combine(other, true); }
    // Each operand may be a row or column vector; preserve the left shape.
    Matrix cross(const Matrix& other) const {
        validate(); other.validate();
        if (!((rows_ == 3 && cols_ == 1) || (rows_ == 1 && cols_ == 3)) ||
            !((other.rows_ == 3 && other.cols_ == 1) || (other.rows_ == 1 && other.cols_ == 3)))
            throw std::invalid_argument("cross product requires three-dimensional vectors");
        const auto& a = values_; const auto& b = other.values_;
        return Matrix(rows_, cols_, {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]});
    }
    Matrix scale(double scalar) const {
        validate(); finite(scalar);
        Matrix result(rows_, cols_);
        for (std::size_t i = 0; i < values_.size(); i++) result.values_[i] = values_[i] * scalar;
        result.validate();
        return result;
    }
    Matrix transpose() const {
        validate();
        Matrix result(cols_, rows_);
        for (std::size_t row = 0; row < rows_; row++)
        for (std::size_t col = 0; col < cols_; col++)
            result.values_[col * rows_ + row] = values_[row * cols_ + col];
        return result;
    }
    Matrix operator*(const Matrix& other) const {
        validate(); other.validate();
        if (cols_ != other.rows_) throw std::invalid_argument("left column count must equal right row count");
        Matrix result(rows_, other.cols_);
        for (std::size_t row = 0; row < rows_; row++)
        for (std::size_t col = 0; col < other.cols_; col++) {
            double sum = 0.0;
            for (std::size_t k = 0; k < cols_; k++) sum += values_[row * cols_ + k] * other.values_[k * other.cols_ + col];
            result.values_[row * other.cols_ + col] = sum;
        }
        result.validate();
        return result;
    }
    double trace() const {
        validate(); square();
        double result = 0.0;
        for (std::size_t i = 0; i < rows_; i++) result += values_[i * cols_ + i];
        return finite(result);
    }
    // Auto never assumes positive definiteness. Explicit Cholesky requires
    // exact symmetry and strictly positive computed pivots.
    double determinant(DeterminantAlgorithm algorithm = DeterminantAlgorithm::Auto) const {
        if (algorithm != DeterminantAlgorithm::Auto && algorithm != DeterminantAlgorithm::Lu &&
            algorithm != DeterminantAlgorithm::Cholesky) throw std::invalid_argument("unknown determinant algorithm");
        validate(); square();
        if (algorithm == DeterminantAlgorithm::Lu) return determinant_lu();
        if (algorithm == DeterminantAlgorithm::Cholesky) return determinant_cholesky();
        if (rows_ == 0) return 1.0;
        if (rows_ == 1) return values_[0];
        bool upper = true, lower = true;
        for (std::size_t row = 0; row < rows_ && (upper || lower); row++)
        for (std::size_t col = 0; col < row; col++) {
            if (values_[row * rows_ + col] != 0.0) upper = false;
            if (values_[col * rows_ + row] != 0.0) lower = false;
        }
        if (upper || lower) {
            DeterminantProduct product;
            for (std::size_t i = 0; i < rows_; i++) product.append(values_[i * rows_ + i]);
            return product.finish();
        }
        double result;
        if (rows_ <= 3 && determinant_small(result)) return result;
        return determinant_lu();
    }
    // Values ascend; corresponding unit eigenvectors are columns of vectors.
    Factorization factor_lu() const;
    Factorization factor_cholesky() const;
    Factorization factor_qr() const;
    struct SpectralDiagnostics { std::size_t rank; double reciprocal_condition, retained_reciprocal_condition; };
    Matrix pseudoinverse(double cutoff=-1) const;
    Matrix solve_ridge(const Matrix& rhs,double lambda) const;
    Matrix solve_minimum_norm(const Matrix& rhs,double cutoff=-1) const;
    SpectralDiagnostics spectral_diagnostics(double cutoff=-1) const;
    struct Svd;
    Svd svd(double tolerance=1e-12, std::size_t max_sweeps=100) const;
    Matrix solve(const Matrix& rhs) const;
    Matrix least_squares(const Matrix& rhs) const;
    GeneralEigenResult eigen_general(std::size_t max_iterations = 1000) const;
    SymmetricEigenResult eigen_symmetric(double tolerance = 1e-12, std::size_t max_sweeps = 50) const;
    std::pair<bool, bool> triangular() const {
        validate();
        bool upper = rows_ == cols_, lower = upper;
        if (upper) for (std::size_t row = 0; row < rows_; row++)
        for (std::size_t col = 0; col < cols_; col++) {
            if (row > col && values_[row * cols_ + col] != 0.0) upper = false;
            if (row < col && values_[row * cols_ + col] != 0.0) lower = false;
        }
        return {upper, lower};
    }
};

struct SymmetricEigenResult {
    std::vector<double> values;
    Matrix vectors;
};

struct GeneralEigenResult {
    std::vector<double> values_real, values_imag;
    Matrix vectors_real, vectors_imag;
};

inline GeneralEigenResult Matrix::eigen_general(std::size_t max_iterations) const {
    validate(); square();
    if (!max_iterations || max_iterations > 100000)
        throw std::invalid_argument("max_iterations must be from 1 to 100000");
    const auto n = rows_;
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("eigendecomposition dimensions overflow");
    GeneralEigenResult result{std::vector<double>(n), std::vector<double>(n), Matrix(n, n), Matrix(n, n)};
    std::vector<double> h(n * n), v(n * n), ort(n), powers(n);
    const int status = la_general_eigen(static_cast<int>(n), values_.data(), max_iterations,
        result.values_real.data(), result.values_imag.data(), result.vectors_real.data(),
        result.vectors_imag.data(), h.data(), v.data(), ort.data(), powers.data());
    if (status == 1) throw EigenConvergenceError();
    if (status == 2) throw std::overflow_error("nonfinite eigendecomposition result");
    return result;
}

inline SymmetricEigenResult Matrix::eigen_symmetric(double tolerance, std::size_t max_sweeps) const {
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || tolerance >= 1.0 || !max_sweeps || max_sweeps > 10000)
        throw std::invalid_argument("expected 0 < tolerance < 1 and 1 <= max_sweeps <= 10000");
    validate(); square();
    const std::size_t n = rows_;
    double maximum = 0.0;
    bool diagonal = true;
    for (std::size_t row = 0; row < n; row++)
    for (std::size_t col = 0; col < n; col++) {
        const double value = values_[row * n + col];
        if (value != values_[col * n + row]) throw NotSymmetricError();
        if (row != col && value != 0.0) diagonal = false;
        maximum = std::max(maximum, std::abs(value));
    }
    SymmetricEigenResult result{std::vector<double>(n), Matrix::identity(n)};
    auto& eigenvectors = result.vectors.values_;
    if (diagonal) {
        for (std::size_t i = 0; i < n; i++) result.values[i] = values_[i * n + i];
    } else {
        auto work = values_;
        int exponent = 0;
        (void)std::frexp(maximum, &exponent);
        double norm = 0.0;
        for (double& value : work) {
            value = std::scalbn(value, -exponent);
            norm = std::hypot(norm, value);
        }
        const double target = tolerance * norm, skip = target / (2.0 * static_cast<double>(n));
        const auto off_norm = [&] {
            double norm = 0.0;
            for (std::size_t row = 0; row < n; row++)
            for (std::size_t col = row + 1; col < n; col++) {
                norm = std::hypot(norm, work[row * n + col]);
                norm = std::hypot(norm, work[row * n + col]);
            }
            return norm;
        };
        std::size_t sweep = 0;
        while (off_norm() > target && sweep < max_sweeps) {
            for (std::size_t p = 0; p < n; p++)
            for (std::size_t q = p + 1; q < n; q++) {
                const double apq = work[p * n + q];
                if (std::abs(apq) <= skip) continue;
                const double app = work[p * n + p], aqq = work[q * n + q];
                const double delta = aqq - app, twice = 2.0 * apq;
                const double t = twice / (delta + std::copysign(std::hypot(delta, twice), delta));
                const double c = 1.0 / std::hypot(1.0, t), s = t * c;
                work[p * n + p] = app - t * apq;
                work[q * n + q] = aqq + t * apq;
                work[p * n + q] = work[q * n + p] = 0.0;
                for (std::size_t k = 0; k < n; k++) {
                    if (k != p && k != q) {
                        const double akp = work[k * n + p], akq = work[k * n + q];
                        work[k * n + p] = work[p * n + k] = c * akp - s * akq;
                        work[k * n + q] = work[q * n + k] = s * akp + c * akq;
                    }
                    const double vkp = eigenvectors[k * n + p], vkq = eigenvectors[k * n + q];
                    eigenvectors[k * n + p] = c * vkp - s * vkq;
                    eigenvectors[k * n + q] = s * vkp + c * vkq;
                }
            }
            sweep++;
        }
        if (off_norm() > target) throw EigenConvergenceError();
        for (std::size_t i = 0; i < n; i++) result.values[i] = finite(std::scalbn(work[i * n + i], exponent));
    }
    for (std::size_t i = 0; i < n; i++) {
        std::size_t smallest = i;
        for (std::size_t j = i + 1; j < n; j++) if (result.values[j] < result.values[smallest]) smallest = j;
        if (smallest != i) {
            std::swap(result.values[i], result.values[smallest]);
            for (std::size_t row = 0; row < n; row++) std::swap(eigenvectors[row * n + i], eigenvectors[row * n + smallest]);
        }
        std::size_t largest = 0;
        for (std::size_t row = 1; row < n; row++)
            if (std::abs(eigenvectors[row * n + i]) > std::abs(eigenvectors[largest * n + i])) largest = row;
        if (eigenvectors[largest * n + i] < 0.0)
            for (std::size_t row = 0; row < n; row++) eigenvectors[row * n + i] = -eigenvectors[row * n + i];
    }
    return result;
}
class Factorization {
    std::unique_ptr<la_factor,decltype(&la_factor_destroy)> factor_{nullptr,la_factor_destroy};
    static void require(int status){
        switch(status){
            case LA_OK:return;
            case LA_MEMORY:throw std::bad_alloc();
            case LA_RANGE:throw std::overflow_error("solver scaling or arithmetic exceeds float64 range");
            case LA_SINGULAR:throw std::invalid_argument("singular matrix: zero computed LU pivot");
            case LA_NOT_SPD:throw std::invalid_argument("Cholesky requires exact symmetry and positive pivots");
            case LA_RANK:throw std::invalid_argument("QR input is numerically rank deficient");
            default:throw std::invalid_argument("invalid factorization or right-hand side shape");
        }
    }
public:
    Factorization(const Matrix& source,int algorithm){
        la_factor *result=nullptr;
        require(la_factor_create(source.rows(),source.cols(),source.values().data(),algorithm,&result));
        factor_.reset(result);
    }
    Matrix solve(const Matrix& rhs) const {
        if(!factor_)throw std::logic_error("factorization was moved");
        if(rhs.rows()!=factor_->rows)throw std::invalid_argument("right-hand side row count must match factorization");
        std::vector<double> values(factor_->cols*rhs.cols());
        require(la_factor_solve(factor_.get(),rhs.rows(),rhs.cols(),rhs.values().data(),values.data(),1));
        return Matrix(factor_->cols,rhs.cols(),values);
    }
    double reciprocal_condition() const {double result=0;require(la_factor_rcond(factor_.get(),&result));return result;}
};
inline Factorization Matrix::factor_lu() const{return Factorization(*this,LA_LU);}
inline Factorization Matrix::factor_cholesky() const{return Factorization(*this,LA_CHOLESKY);}
inline Factorization Matrix::factor_qr() const{return Factorization(*this,LA_QR);}
inline Matrix Matrix::solve(const Matrix& rhs) const{return factor_lu().solve(rhs);}
inline Matrix Matrix::least_squares(const Matrix& rhs) const{return factor_qr().solve(rhs);}

struct Matrix::Svd { Matrix u; std::vector<double> values; Matrix vt; };
inline Matrix::Svd Matrix::svd(double tolerance,std::size_t max_sweeps) const {
    auto k=std::min(rows_,cols_);
    Svd result{Matrix(rows_,k),std::vector<double>(k),Matrix(k,cols_)};
    int status=la_svd(rows_,cols_,values_.data(),tolerance,max_sweeps,result.u.data(),result.values.data(),result.vt.data());
    if(status==1)throw std::invalid_argument("invalid SVD options");
    if(status==2)throw std::bad_alloc();
    if(status==3)throw std::overflow_error("SVD outside float64 range");
    if(status==4)throw EigenConvergenceError();
    return result;
}
inline double inverse_cutoff(std::size_t m,std::size_t n,double cutoff) {
    if(cutoff==-1)cutoff=static_cast<double>(std::max(m,n))*std::numeric_limits<double>::epsilon();
    if(!std::isfinite(cutoff)||cutoff<0||cutoff>1)throw std::invalid_argument("invalid relative cutoff");
    return cutoff;
}
inline Matrix Matrix::pseudoinverse(double cutoff) const {
    cutoff=inverse_cutoff(rows_,cols_,cutoff);auto r=svd();Matrix out(cols_,rows_);
    if(la_apply_inverse(rows_,cols_,r.values.size(),r.u.data(),r.values.data(),r.vt.data(),cutoff,nullptr,rows_,0,out.data(),0))throw std::overflow_error("SVD inverse outside float64 range");
    return out;
}
inline Matrix Matrix::solve_minimum_norm(const Matrix& rhs,double cutoff) const {
    rhs.validate();
    if(rhs.rows_!=rows_)throw std::invalid_argument("incompatible right-hand side");
    cutoff=inverse_cutoff(rows_,cols_,cutoff);auto r=svd();Matrix out(cols_,rhs.cols_);
    if(la_apply_inverse(rows_,cols_,r.values.size(),r.u.data(),r.values.data(),r.vt.data(),cutoff,rhs.values_.data(),rhs.cols_,1,out.data(),0))throw std::overflow_error("SVD solve outside float64 range");
    return out;
}
inline Matrix Matrix::solve_ridge(const Matrix& rhs,double lambda) const {
    if(!std::isfinite(lambda)||lambda<0)throw std::invalid_argument("lambda must be finite and nonnegative");
    if(lambda==0)return solve_minimum_norm(rhs);
    rhs.validate();if(rhs.rows_!=rows_)throw std::invalid_argument("incompatible right-hand side");
    auto r=svd();Matrix out(cols_,rhs.cols_);
    if(la_apply_inverse(rows_,cols_,r.values.size(),r.u.data(),r.values.data(),r.vt.data(),0,rhs.values_.data(),rhs.cols_,1,out.data(),lambda))throw std::overflow_error("Ridge solve outside float64 range");
    return out;
}
inline Matrix::SpectralDiagnostics Matrix::spectral_diagnostics(double cutoff) const {
    cutoff=inverse_cutoff(rows_,cols_,cutoff);auto r=svd();auto k=r.values.size();auto rank=la_spectral_rank(k,r.values.data(),cutoff);
    return {rank,k&&r.values[0]?r.values[k-1]/r.values[0]:0,rank?r.values[rank-1]/r.values[0]:0};
}
}
#endif
