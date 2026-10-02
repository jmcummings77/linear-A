#include "matrix.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using linear_a::Matrix;
using linear_a::DeterminantAlgorithm;
#define CHECK(condition) do { if (!(condition)) throw std::runtime_error("test failed: " #condition); } while (false)
template<class Exception, class Operation> void throws(Operation operation) {
    try { operation(); } catch (const Exception&) { return; }
    throw std::runtime_error("expected exception was not thrown");
}
static std::int64_t exact_determinant(const std::vector<double>& values, std::size_t n) {
    if (!n) return 1;
    std::int64_t result = 0;
    for (std::size_t col = 0; col < n; col++) {
        std::vector<double> minor;
        for (std::size_t row = 1; row < n; row++)
        for (std::size_t j = 0; j < n; j++) if (j != col) minor.push_back(values[row * n + j]);
        auto term = static_cast<std::int64_t>(values[col]) * exact_determinant(minor, n - 1);
        result += col % 2 ? -term : term;
    }
    return result;
}
static void determinant_tests() {
    const auto algorithms = {DeterminantAlgorithm::Auto, DeterminantAlgorithm::Lu, DeterminantAlgorithm::Cholesky};
    std::uint32_t state = 17;
    for (std::size_t n = 0; n <= 4; n++)
    for (std::size_t sample = 0; sample < 80; sample++) {
        std::vector<double> values(n * n);
        for (double& value : values) {
            state = state * 1664525u + 1013904223u;
            value = static_cast<double>(state % 7) - 3;
        }
        Matrix matrix(n, n, values);
        const double expected = static_cast<double>(exact_determinant(values, n));
        for (auto algorithm : {DeterminantAlgorithm::Auto, DeterminantAlgorithm::Lu}) {
            CHECK(std::abs(matrix.determinant(algorithm) - expected) <= 1e-10 * std::max(1.0, std::abs(expected)));
            CHECK(matrix.values() == values);
        }
    }
    for (std::size_t sample = 0; sample < 40; sample++) {
        std::vector<double> lower(16, 0), values(16, 0);
        double expected = 1;
        for (std::size_t row = 0; row < 4; row++)
        for (std::size_t col = 0; col <= row; col++) {
            state = state * 1664525u + 1013904223u;
            lower[row * 4 + col] = row == col ? static_cast<double>(state % 4 + 1) : static_cast<double>(state % 5) - 2;
            if (row == col) expected *= lower[row * 4 + col] * lower[row * 4 + col];
        }
        for (std::size_t row = 0; row < 4; row++)
        for (std::size_t col = 0; col < 4; col++)
        for (std::size_t k = 0; k < 4; k++) values[row * 4 + col] += lower[row * 4 + k] * lower[col * 4 + k];
        const Matrix matrix(4, 4, values);
        for (auto algorithm : algorithms) {
            CHECK(std::abs(matrix.determinant(algorithm) - expected) < 1e-10 * expected);
            CHECK(matrix.values() == values);
        }
    }
    struct Example { std::size_t n; std::vector<double> values; double expected; };
    const double tiny = std::numeric_limits<double>::denorm_min(), epsilon = std::numeric_limits<double>::epsilon();
    const Example examples[] = {
        {2, {1e200, 1e200, 1e200, 1e200}, 0},
        {2, {1e-200, 1e-200, 3e-124, 6e-124}, tiny},
        {3, {1e-200, 1e-200, 0, 3e-124, 6e-124, 0, 0, 0, 1}, tiny},
        {2, {1, 1, 1, 1 + epsilon}, epsilon},
        {3, {1e200, 1e200, 0, 1e200, 1e200, 0, 0, 0, 1}, 0},
    };
    for (const auto& example : examples) {
        const Matrix matrix(example.n, example.n, example.values);
        for (auto algorithm : {DeterminantAlgorithm::Auto, DeterminantAlgorithm::Lu})
            CHECK(matrix.determinant(algorithm) == example.expected);
    }
    for (const auto& values : std::vector<std::vector<double>>{{1e308, 1e308, 1e-308, 2e-308}, {1e308, 1e-308, 1e308, 2e-308}}) {
        const Matrix matrix(2, 2, values);
        for (auto algorithm : {DeterminantAlgorithm::Auto, DeterminantAlgorithm::Lu})
            CHECK(std::abs(matrix.determinant(algorithm) - 1) < 1e-14);
    }
    for (const auto& exponents : std::vector<std::vector<int>>{{800, 800, -800, -800}, {-800, -800, 800, 800}, {800, -800, 800, -800}}) {
        Matrix diagonal(4, 4);
        for (std::size_t i = 0; i < 4; i++) diagonal.set(i, i, std::scalbn(1.0, exponents[i]));
        for (auto algorithm : algorithms) CHECK(diagonal.determinant(algorithm) == 1);
        diagonal.set(3, 3, 0);
        CHECK(diagonal.determinant() == 0);
    }
    const double maximum = std::numeric_limits<double>::max();
    for (auto algorithm : algorithms) {
        CHECK(Matrix().determinant(algorithm) == 1);
        CHECK(Matrix(1, 1, {tiny}).determinant(algorithm) == tiny);
        throws<std::overflow_error>([&] { (void)Matrix(2, 2, {maximum, 0, 0, 2}).determinant(algorithm); });
        throws<std::invalid_argument>([&] { (void)Matrix(2, 3).determinant(algorithm); });
        Matrix corrupt(1, 1); corrupt.data()[0] = std::numeric_limits<double>::infinity();
        throws<std::overflow_error>([&] { (void)corrupt.determinant(algorithm); });
    }
    CHECK(std::abs(Matrix(2, 2, {4e300, 2, 2, 3e-300}).determinant(DeterminantAlgorithm::Cholesky) - 8) < 1e-13);
    for (const auto& values : std::vector<std::vector<double>>{{2, 100, 1, 2}, {1, 2, 2, 1}, {1, 1, 1, 1}, {-1, 0, 0, -1}, {2, 1, 1 + epsilon, 2}}) {
        const Matrix matrix(2, 2, values);
        throws<std::domain_error>([&] { (void)matrix.determinant(DeterminantAlgorithm::Cholesky); });
        CHECK(matrix.values() == values);
    }
    throws<std::invalid_argument>([] { (void)Matrix().determinant(static_cast<DeterminantAlgorithm>(99)); });
}

static void eigen_invariants(const Matrix& a, const linear_a::SymmetricEigenResult& result) {
    const auto n = a.rows();
    CHECK(result.values.size() == n && result.vectors.rows() == n && result.vectors.cols() == n);
    double scale = 0;
    for (double value : a.values()) scale = std::max(scale, std::abs(value));
    if (scale == 0) scale = 1;
    for (std::size_t col = 0; col < n; col++) {
        CHECK(std::isfinite(result.values[col]));
        if (col) CHECK(result.values[col - 1] <= result.values[col]);
        std::size_t largest = 0;
        for (std::size_t row = 0; row < n; row++) {
            double av = 0, dot = 0;
            for (std::size_t k = 0; k < n; k++) {
                av += (a.get(row, k) / scale) * result.vectors.get(k, col);
                dot += result.vectors.get(k, row) * result.vectors.get(k, col);
            }
            CHECK(std::abs(av - result.vectors.get(row, col) * (result.values[col] / scale)) <= 2e-11 * n);
            CHECK(std::abs(dot - (row == col ? 1.0 : 0.0)) <= 2e-13 * n);
            if (std::abs(result.vectors.get(row, col)) > std::abs(result.vectors.get(largest, col))) largest = row;
        }
        CHECK(result.vectors.get(largest, col) >= 0);
    }
}
static void eigen_tests() {
    for (std::size_t n = 0; n <= 12; n++) {
        Matrix a(n, n);
        for (std::size_t row = 0; row < n; row++) {
            a.set(row, row, 2.25);
            if (row + 1 < n) { a.set(row, row + 1, -1); a.set(row + 1, row, -1); }
        }
        const auto before = a.values();
        auto result = a.eigen_symmetric();
        eigen_invariants(a, result);
        for (std::size_t i = 0; i < n; i++) CHECK(std::abs(result.values[i] - (2.25 - 2 * std::cos((i + 1) * std::acos(-1.0) / (n + 1)))) < 2e-12);
        CHECK(before == a.values());
    }
    const double tiny = std::numeric_limits<double>::denorm_min(), maximum = std::numeric_limits<double>::max();
    for (double scale : {1.0, 1e300, 1e-300, tiny}) {
        const Matrix a(2, 2, {2 * scale, scale, scale, 2 * scale});
        auto result = a.eigen_symmetric(); eigen_invariants(a, result);
        CHECK(std::abs(result.values[0] / scale - 1) < 1e-14 && std::abs(result.values[1] / scale - 3) < 1e-14);
    }
    const Matrix diagonal(3, 3, {maximum, 0, 0, 0, tiny, 0, 0, 0, -maximum});
    auto result = diagonal.eigen_symmetric(); eigen_invariants(diagonal, result);
    CHECK(result.values == std::vector<double>({-maximum, tiny, maximum}));
    eigen_invariants(Matrix::identity(4), Matrix::identity(4).eigen_symmetric());
    const Matrix dense(4, 4, {4,1,2,3, 1,5,1,2, 2,1,6,1, 3,2,1,7});
    throws<linear_a::EigenConvergenceError>([&] { (void)dense.eigen_symmetric(1e-15, 1); });
    for (double tolerance : {0.0, -1.0, 1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        throws<std::invalid_argument>([&] { (void)dense.eigen_symmetric(tolerance); });
    throws<std::invalid_argument>([&] { (void)dense.eigen_symmetric(1e-12, 0); });
    throws<std::invalid_argument>([&] { (void)dense.eigen_symmetric(1e-12, 10001); });
    throws<linear_a::NotSymmetricError>([] { (void)Matrix(2, 2, {1, 2, 2 + 1e-15, 1}).eigen_symmetric(); });
    throws<std::invalid_argument>([] { (void)Matrix(2, 3).eigen_symmetric(); });
    throws<std::overflow_error>([&] { (void)Matrix(2, 2, {maximum, maximum, maximum, maximum}).eigen_symmetric(); });
    Matrix corrupt(1, 1); corrupt.data()[0] = std::numeric_limits<double>::quiet_NaN();
    throws<std::overflow_error>([&] { (void)corrupt.eigen_symmetric(); });
    std::uint32_t state = 31;
    for (std::size_t sample = 0; sample < 30; sample++) {
        Matrix a(7, 7);
        for (std::size_t row = 0; row < 7; row++) for (std::size_t col = row; col < 7; col++) {
            state = state * 1664525u + 1013904223u;
            const double value = static_cast<double>(state % 101) - 50;
            a.set(row, col, value); a.set(col, row, value);
        }
        eigen_invariants(a, a.eigen_symmetric());
    }
}


static void rotation_invariants(const Matrix& r) {
    const auto identity = r.transpose() * r;
    for (std::size_t row = 0; row < r.rows(); row++) for (std::size_t col = 0; col < r.cols(); col++)
        CHECK(std::abs(identity.get(row,col)-(row == col ? 1.0 : 0.0)) < 2e-14);
    CHECK(std::abs(r.determinant()-1) < 2e-14);
}
static void vector_rotation_tests() {
    for (std::size_t left : {1,3}) for (std::size_t right : {1,3}) {
        const Matrix a(left,3/left,{2,3,4}), b(right,3/right,{5,6,7});
        auto result = a.cross(b), reverse = b.cross(a);
        CHECK(result.rows() == left && result.cols() == 3/left);
        CHECK(result.values() == std::vector<double>({-3,6,-3}));
        CHECK(reverse.values() == std::vector<double>({3,-6,3}));
        CHECK(a.cross(a).values() == std::vector<double>({0,0,0}));
        result.set(0,0,22); CHECK(a.values() == std::vector<double>({2,3,4}));
    }
    const Matrix ex(3,1,{1,0,0}), ey(3,1,{0,1,0});
    CHECK(ex.cross(ey).values() == std::vector<double>({0,0,1}));
    throws<std::invalid_argument>([&] { (void)Matrix(3,3).cross(ex); });
    throws<std::invalid_argument>([&] { (void)ex.cross(Matrix(1,2)); });
    const double maximum = std::numeric_limits<double>::max(), infinity = std::numeric_limits<double>::infinity();
    throws<std::overflow_error>([&] { (void)Matrix(3,1,{maximum,0,0}).cross(Matrix(3,1,{0,2,0})); });
    Matrix corrupt(3,1); corrupt.data()[0] = infinity;
    throws<std::overflow_error>([&] { (void)ex.cross(corrupt); });
    const double quarter = std::acos(-1.0)/2;
    CHECK(std::abs((Matrix::rotation_2d(quarter)*Matrix(2,1,{1,0})).get(1,0)-1) < 1e-15);
    CHECK(std::abs((Matrix::rotation_x(quarter)*ey).get(2,0)-1) < 1e-15);
    CHECK(std::abs((Matrix::rotation_y(quarter)*ex).get(2,0)+1) < 1e-15);
    CHECK(std::abs((Matrix::rotation_z(quarter)*ex).get(1,0)-1) < 1e-15);
    for (double angle : {0.0,quarter,-std::acos(-1.0),4*std::acos(-1.0),1e-8,.5,maximum}) {
        rotation_invariants(Matrix::rotation_2d(angle));
        const Matrix fixed[] = {Matrix::rotation_x(angle),Matrix::rotation_y(angle),Matrix::rotation_z(angle)};
        for (std::size_t index = 0; index < 3; index++) {
            Matrix axis(3,1); axis.set(index,0,1);
            auto general = Matrix::rotation_axis_angle(axis,angle);
            rotation_invariants(fixed[index]); rotation_invariants(general);
            for (std::size_t i = 0; i < 9; i++) CHECK(std::abs(general.values()[i]-fixed[index].values()[i]) < 2e-15);
        }
    }
    auto normal = Matrix::rotation_axis_angle(Matrix(1,3,{1,1,1}),.5);
    for (double scale : {1.0,1e-300,1e300,std::numeric_limits<double>::denorm_min()}) {
        const Matrix axis(3,1,{scale,scale,scale});
        auto r = Matrix::rotation_axis_angle(axis,.5); rotation_invariants(r);
        for (std::size_t i = 0; i < 9; i++) CHECK(std::abs(r.values()[i]-normal.values()[i]) < 2e-15);
        CHECK(axis.values() == std::vector<double>({scale,scale,scale}));
    }
    const Matrix axis(3,1,{1,2,3}), a(3,1,{2,3,4}), b(3,1,{5,6,7});
    const auto r = Matrix::rotation_axis_angle(axis,.5);
    auto axis_result = r*axis, cross_rotated = (r*a).cross(r*b), rotated_cross = r*a.cross(b);
    for (std::size_t i = 0; i < 3; i++) {
        CHECK(std::abs(axis_result.values()[i]-axis.values()[i]) < 2e-14);
        CHECK(std::abs(cross_rotated.values()[i]-rotated_cross.values()[i]) < 2e-14);
    }
    const auto inverse = Matrix::rotation_axis_angle(axis.scale(-1),.5);
    const auto transposed = r.transpose();
    for (std::size_t i = 0; i < 9; i++) CHECK(std::abs(inverse.values()[i]-transposed.values()[i]) < 2e-15);
    CHECK(std::abs(Matrix::rotation_axis_angle(Matrix(3,1,{1,1,0}),1e-8).get(0,1)/2.5e-17-1) < 1e-14);
    throws<std::invalid_argument>([] { (void)Matrix::rotation_axis_angle(Matrix(3,1),.5); });
    throws<std::invalid_argument>([] { (void)Matrix::rotation_axis_angle(Matrix(3,3),.5); });
    throws<std::overflow_error>([&] { (void)Matrix::rotation_axis_angle(corrupt,.5); });
    throws<std::overflow_error>([&] { (void)Matrix::rotation_axis_angle(axis,infinity); });
    throws<std::overflow_error>([&] { (void)Matrix::rotation_2d(infinity); });
    throws<std::overflow_error>([&] { (void)Matrix::rotation_x(infinity); });
    throws<std::overflow_error>([&] { (void)Matrix::rotation_y(infinity); });
    throws<std::overflow_error>([&] { (void)Matrix::rotation_z(infinity); });
}


static void general_tests() {
    const auto a = Matrix(2,2,{0,-1,1,0});
    const auto pair = a.eigen_general();
    CHECK(pair.values_real == std::vector<double>({0,0}));
    CHECK(pair.values_imag == std::vector<double>({-1,1}));
    CHECK(a.values() == std::vector<double>({0,-1,1,0}));
    for (std::size_t col=0; col<2; col++) {
        double norm=0;
        for (std::size_t row=0; row<2; row++) {
            double real=0,imag=0;
            for (std::size_t k=0; k<2; k++) {
                real+=a.get(row,k)*pair.vectors_real.get(k,col);
                imag+=a.get(row,k)*pair.vectors_imag.get(k,col);
            }
            CHECK(std::abs(real+pair.values_imag[col]*pair.vectors_imag.get(row,col))<1e-13);
            CHECK(std::abs(imag-pair.values_imag[col]*pair.vectors_real.get(row,col))<1e-13);
            norm+=std::pow(pair.vectors_real.get(row,col),2)+std::pow(pair.vectors_imag.get(row,col),2);
        }
        CHECK(std::abs(norm-1)<1e-13);
    }
    auto copy=pair; copy.vectors_real.set(0,0,9);
    CHECK(pair.vectors_real.get(0,0)!=9);
    const auto empty=Matrix().eigen_general(); CHECK(empty.values_real.empty() && empty.vectors_real.rows()==0);
    auto jordan=Matrix(3,3,{1,1,0,0,1,1,0,0,1}).eigen_general();
    CHECK(jordan.values_real==std::vector<double>({1,1,1}));
    auto mixed=Matrix(2,2,{0,1e300,-1e-300,0}).eigen_general();
    CHECK(std::abs(mixed.values_imag[0]+1)<1e-14 && std::abs(mixed.values_imag[1]-1)<1e-14);
    throws<std::invalid_argument>([&] { (void)a.eigen_general(0); });
    throws<std::invalid_argument>([&] { (void)a.eigen_general(100001); });
    throws<std::invalid_argument>([] { (void)Matrix(2,3).eigen_general(); });
    throws<linear_a::EigenConvergenceError>([] { (void)Matrix(3,3,{1,2,3,4,5,6,7,8,10}).eigen_general(1); });
    const auto maximum=std::numeric_limits<double>::max();
    throws<std::overflow_error>([&] { (void)Matrix(2,2,{maximum,maximum,maximum,maximum}).eigen_general(); });
    auto corrupt=Matrix::identity(2); corrupt.data()[0]=std::numeric_limits<double>::quiet_NaN();
    throws<std::overflow_error>([&] { (void)corrupt.eigen_general(); });
}

int main() {
    general_tests();
    determinant_tests();
    eigen_tests();
    vector_rotation_tests();
    Matrix a(2, 3, {1, 2, 3, 4, 5, 6});
    auto copy = a;
    copy.set(0, 0, 99);
    CHECK(a.get(0, 0) == 1);
    auto moved = std::move(copy);
    CHECK(moved.get(0, 0) == 99 && copy.rows() == 0 && copy.cols() == 0);
    CHECK(a.transpose().values() == std::vector<double>({1, 4, 2, 5, 3, 6}));
    CHECK(a.transpose().rows() == 3 && a.transpose().cols() == 2);
    CHECK(a.row(1) == std::vector<double>({4, 5, 6}));
    CHECK(a.column(1) == std::vector<double>({2, 5}));
    auto row = a.row(0); row[0] = 100;
    CHECK(a.get(0, 0) == 1);
    CHECK((a + a).values() == std::vector<double>({2, 4, 6, 8, 10, 12}));
    CHECK((a - a).values() == std::vector<double>(6, 0));
    CHECK(a.scale(0.5).values() == std::vector<double>({0.5, 1, 1.5, 2, 2.5, 3}));
    Matrix b(3, 2, {7, 8, 9, 10, 11, 12});
    CHECK((a * b).values() == std::vector<double>({58, 64, 139, 154}));
    CHECK(a.get(0, 0) == 1 && b.get(2, 1) == 12);
    CHECK(!a.triangular().first && !a.triangular().second);
    throws<std::invalid_argument>([&] { (void)(a + b); });
    throws<std::invalid_argument>([&] { (void)(a * a); });
    throws<std::invalid_argument>([&] { (void)a.trace(); });
    throws<std::invalid_argument>([&] { (void)a.determinant(); });
    throws<std::out_of_range>([&] { (void)a.get(2, 0); });
    throws<std::out_of_range>([&] { (void)a.row(2); });
    throws<std::out_of_range>([&] { (void)a.column(3); });
    throws<std::invalid_argument>([] { Matrix(2, 2, {1}); });
    throws<std::length_error>([] { Matrix(std::numeric_limits<std::size_t>::max(), 2); });
    throws<std::overflow_error>([&] { a.set(0, 0, std::numeric_limits<double>::infinity()); });
    CHECK(a.get(0, 0) == 1);
    throws<std::overflow_error>([] { Matrix(1, 1, {1e308}).scale(2); });
    Matrix square(3, 3, {6, 1, 1, 4, -2, 5, 2, 8, 7});
    CHECK(std::abs(square.determinant() + 306) < 1e-10);
    CHECK(square.trace() == 11 && square.get(0, 0) == 6);
    CHECK(Matrix(2, 2, {0, 1, 1, 0}).determinant() == -1);
    CHECK(Matrix(2, 2, {1, 2, 2, 4}).determinant() == 0);
    CHECK(Matrix(1, 1, {-7}).determinant() == -7);
    CHECK(Matrix().determinant() == 1 && Matrix().trace() == 0);
    CHECK((Matrix(2, 0) * Matrix(0, 3)).values() == std::vector<double>(6, 0));
    auto identity = Matrix::identity(3);
    CHECK(identity.determinant() == 1 && identity.triangular().first && identity.triangular().second);
    identity.set(0, 1, 3);
    CHECK(identity.triangular().first && !identity.triangular().second);
    std::cout << "C++ matrix tests passed\n";
}
