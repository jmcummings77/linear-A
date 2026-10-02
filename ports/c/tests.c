#include "matrix.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { fprintf(stderr, "test failed at line %d: %s\n", __LINE__, #condition); return 1; } } while (0)
#define OK(operation) CHECK((operation) == M_OK)

/* Small integer expansion is independent of all production algorithms. */
static int64_t exact_determinant(const double *values, size_t n) {
    if (!n) return 1;
    int64_t result = 0;
    for (size_t col = 0; col < n; col++) {
        double minor[9];
        size_t next = 0;
        for (size_t row = 1; row < n; row++)
        for (size_t j = 0; j < n; j++)
            if (j != col) minor[next++] = values[row * n + j];
        int64_t term = (int64_t)values[col] * exact_determinant(minor, n - 1);
        result += col % 2 ? -term : term;
    }
    return result;
}

static int determinant_tests(void) {
    matrix a = {0};
    double value;
    uint32_t state = 17;
    for (size_t n = 0; n <= 4; n++)
    for (size_t sample = 0; sample < 80; sample++) {
        double values[16];
        for (size_t i = 0; i < n * n; i++) {
            state = state * 1664525u + 1013904223u;
            values[i] = (double)(state % 7) - 3;
        }
        OK(m_from_array(n, n, values, &a));
        double expected = (double)exact_determinant(values, n);
        for (int algorithm = M_DETERMINANT_AUTO; algorithm <= M_DETERMINANT_LU; algorithm++) {
            OK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)algorithm, &value));
            CHECK(fabs(value - expected) <= 1e-10 * fmax(1.0, fabs(expected)));
            CHECK(!n || memcmp(values, a.values, n * n * sizeof(double)) == 0);
        }
        m_free(&a);
    }
    /* Independently known det(L L^T) = product(diag(L))^2. */
    for (size_t sample = 0; sample < 40; sample++) {
        double lower[16] = {0}, values[16] = {0}, expected = 1;
        for (size_t row = 0; row < 4; row++)
        for (size_t col = 0; col <= row; col++) {
            state = state * 1664525u + 1013904223u;
            lower[row * 4 + col] = row == col ? (double)(state % 4 + 1) : (double)(state % 5) - 2;
            if (row == col) expected *= lower[row * 4 + col] * lower[row * 4 + col];
        }
        for (size_t row = 0; row < 4; row++)
        for (size_t col = 0; col < 4; col++)
        for (size_t k = 0; k < 4; k++) values[row * 4 + col] += lower[row * 4 + k] * lower[col * 4 + k];
        OK(m_from_array(4, 4, values, &a));
        for (int algorithm = M_DETERMINANT_AUTO; algorithm <= M_DETERMINANT_CHOLESKY; algorithm++) {
            OK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)algorithm, &value));
            CHECK(fabs(value - expected) <= 1e-10 * expected);
            CHECK(memcmp(values, a.values, sizeof(values)) == 0);
        }
        m_free(&a);
    }
    const double extreme[][9] = {
        {1e200, 1e200, 1e200, 1e200},
        {1e308, 1e308, 1e-308, 2e-308},
        {1e308, 1e-308, 1e308, 2e-308},
        {1e-200, 1e-200, 3e-124, 6e-124},
        {1e-200, 1e-200, 0, 3e-124, 6e-124, 0, 0, 0, 1},
        {1, 1, 1, 1 + DBL_EPSILON},
        {1e200, 1e200, 0, 1e200, 1e200, 0, 0, 0, 1},
    };
    const size_t sizes[] = {2, 2, 2, 2, 3, 2, 3};
    const double expected[] = {0, 1, 1, DBL_TRUE_MIN, DBL_TRUE_MIN, DBL_EPSILON, 0};
    for (size_t i = 0; i < sizeof(sizes) / sizeof(*sizes); i++) {
        OK(m_from_array(sizes[i], sizes[i], extreme[i], &a));
        for (int algorithm = M_DETERMINANT_AUTO; algorithm <= M_DETERMINANT_LU; algorithm++) {
            OK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)algorithm, &value));
            if (i == 1 || i == 2) CHECK(fabs(value - expected[i]) < 1e-14);
            else CHECK(value == expected[i]);
        }
        m_free(&a);
    }
    /* Powers of two let us demand an exact product in several pivot orders. */
    const int orders[][4] = {{800, 800, -800, -800}, {-800, -800, 800, 800}, {800, -800, 800, -800}};
    for (size_t order = 0; order < 3; order++) {
        OK(m_create(4, 4, &a));
        for (size_t i = 0; i < 4; i++) a.values[i * 4 + i] = scalbn(1.0, orders[order][i]);
        for (int algorithm = M_DETERMINANT_AUTO; algorithm <= M_DETERMINANT_CHOLESKY; algorithm++) {
            OK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)algorithm, &value));
            CHECK(value == 1.0);
        }
        a.values[3 * 4 + 3] = 0.0;
        OK(m_determinant(&a, &value)); CHECK(value == 0.0);
        m_free(&a);
    }
    const double mixed_spd[] = {4e300, 2, 2, 3e-300};
    OK(m_from_array(2, 2, mixed_spd, &a));
    OK(m_determinant_with_algorithm(&a, M_DETERMINANT_CHOLESKY, &value));
    CHECK(fabs(value - 8) < 1e-13);
    m_free(&a);
    const double not_spd[][4] = {{2, 100, 1, 2}, {1, 2, 2, 1}, {1, 1, 1, 1}, {-1, 0, 0, -1}, {2, 1, 1 + DBL_EPSILON, 2}};
    for (size_t i = 0; i < sizeof(not_spd) / sizeof(*not_spd); i++) {
        OK(m_from_array(2, 2, not_spd[i], &a));
        value = 123;
        CHECK(m_determinant_with_algorithm(&a, M_DETERMINANT_CHOLESKY, &value) == M_NOT_POSITIVE_DEFINITE);
        CHECK(value == 123 && memcmp(not_spd[i], a.values, sizeof(not_spd[i])) == 0);
        m_free(&a);
    }
    OK(m_create(0, 0, &a));
    OK(m_determinant_with_algorithm(&a, M_DETERMINANT_CHOLESKY, &value)); CHECK(value == 1);
    value = 123;
    CHECK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)99, &value) == M_ALGORITHM && value == 123);
    CHECK(m_determinant_with_algorithm(&a, M_DETERMINANT_AUTO, NULL) == M_ARGUMENT);
    m_free(&a);
    const double tiny[] = {DBL_TRUE_MIN};
    OK(m_from_array(1, 1, tiny, &a));
    for (int algorithm = M_DETERMINANT_AUTO; algorithm <= M_DETERMINANT_CHOLESKY; algorithm++) {
        OK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)algorithm, &value));
        CHECK(value == DBL_TRUE_MIN);
    }
    a.values[0] = INFINITY;
    for (int algorithm = M_DETERMINANT_AUTO; algorithm <= M_DETERMINANT_CHOLESKY; algorithm++) {
        value = 123;
        CHECK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)algorithm, &value) == M_NONFINITE && value == 123);
    }
    m_free(&a);
    const double overflow[] = {DBL_MAX, 0, 0, 2};
    OK(m_from_array(2, 2, overflow, &a));
    for (int algorithm = M_DETERMINANT_AUTO; algorithm <= M_DETERMINANT_CHOLESKY; algorithm++) {
        value = 123;
        CHECK(m_determinant_with_algorithm(&a, (matrix_determinant_algorithm)algorithm, &value) == M_NONFINITE && value == 123);
    }
    m_free(&a);
    return 0;
}


static int eigen_invariants(const matrix *a, const matrix *values, const matrix *vectors) {
    const size_t n = a->rows;
    CHECK(values->rows == n && values->cols == 1 && vectors->rows == n && vectors->cols == n);
    double scale = 0.0;
    for (size_t i = 0; i < n * n; i++) scale = fmax(scale, fabs(a->values[i]));
    if (scale == 0.0) scale = 1.0;
    for (size_t col = 0; col < n; col++) {
        CHECK(isfinite(values->values[col]));
        if (col) CHECK(values->values[col - 1] <= values->values[col]);
        size_t largest = 0;
        for (size_t row = 0; row < n; row++) {
            double av = 0.0, vd = vectors->values[row * n + col] * (values->values[col] / scale);
            for (size_t k = 0; k < n; k++) av += (a->values[row * n + k] / scale) * vectors->values[k * n + col];
            CHECK(fabs(av - vd) <= 2e-11 * (double)n);
            if (fabs(vectors->values[row * n + col]) > fabs(vectors->values[largest * n + col])) largest = row;
            double dot = 0.0;
            for (size_t k = 0; k < n; k++) dot += vectors->values[k * n + row] * vectors->values[k * n + col];
            CHECK(fabs(dot - (row == col ? 1.0 : 0.0)) <= 2e-13 * (double)n);
        }
        CHECK(vectors->values[largest * n + col] >= 0.0);
    }
    return 0;
}

static int eigen_tests(void) {
    matrix a = {0}, values = {0}, vectors = {0};
    for (size_t n = 0; n <= 12; n++) {
        OK(m_create(n, n, &a));
        for (size_t row = 0; row < n; row++) {
            a.values[row * n + row] = 2.25;
            if (row + 1 < n) a.values[row * n + row + 1] = a.values[(row + 1) * n + row] = -1;
        }
        matrix copy = {0}; OK(m_copy(&a, &copy));
        OK(m_eigen_symmetric(&a, &values, &vectors));
        CHECK(eigen_invariants(&a, &values, &vectors) == 0);
        for (size_t i = 0; i < n; i++) CHECK(fabs(values.values[i] - (2.25 - 2 * cos((double)(i + 1) * acos(-1.0) / (double)(n + 1)))) < 2e-12);
        CHECK(!n || memcmp(a.values, copy.values, n * n * sizeof(double)) == 0);
        m_free(&copy); m_free(&a); m_free(&values); m_free(&vectors);
    }
    const double scales[] = {1.0, 1e300, 1e-300, DBL_TRUE_MIN};
    for (size_t i = 0; i < 4; i++) {
        const double input[] = {2 * scales[i], scales[i], scales[i], 2 * scales[i]};
        OK(m_from_array(2, 2, input, &a)); OK(m_eigen_symmetric(&a, &values, &vectors));
        CHECK(eigen_invariants(&a, &values, &vectors) == 0);
        CHECK(fabs(values.values[0] / scales[i] - 1.0) < 1e-14);
        CHECK(fabs(values.values[1] / scales[i] - 3.0) < 1e-14);
        m_free(&a); m_free(&values); m_free(&vectors);
    }
    uint32_t state = 31;
    for (size_t sample = 0; sample < 30; sample++) {
        OK(m_create(7, 7, &a));
        for (size_t row = 0; row < 7; row++) for (size_t col = row; col < 7; col++) {
            state = state * 1664525u + 1013904223u;
            a.values[row * 7 + col] = a.values[col * 7 + row] = (double)(state % 101) - 50;
        }
        OK(m_eigen_symmetric(&a, &values, &vectors)); CHECK(eigen_invariants(&a, &values, &vectors) == 0);
        m_free(&a); m_free(&values); m_free(&vectors);
    }
    const double diagonal[] = {DBL_MAX, 0, 0, 0, DBL_TRUE_MIN, 0, 0, 0, -DBL_MAX};
    OK(m_from_array(3, 3, diagonal, &a)); OK(m_eigen_symmetric(&a, &values, &vectors));
    CHECK(values.values[0] == -DBL_MAX && values.values[1] == DBL_TRUE_MIN && values.values[2] == DBL_MAX);
    CHECK(eigen_invariants(&a, &values, &vectors) == 0);
    CHECK(m_eigen_symmetric(&a, &values, &vectors) == M_ARGUMENT);
    m_free(&a); m_free(&values); m_free(&vectors);
    OK(m_identity(4, &a)); OK(m_eigen_symmetric(&a, &values, &vectors));
    CHECK(eigen_invariants(&a, &values, &vectors) == 0); m_free(&a); m_free(&values); m_free(&vectors);
    const double dense[] = {4,1,2,3, 1,5,1,2, 2,1,6,1, 3,2,1,7};
    OK(m_from_array(4, 4, dense, &a));
    CHECK(m_eigen_symmetric_with_options(&a, 1e-15, 1, &values, &vectors) == M_NO_CONVERGENCE);
    CHECK(values.values == NULL && vectors.values == NULL && values.cols == 0 && vectors.cols == 0);
    CHECK(memcmp(a.values, dense, sizeof(dense)) == 0);
    const double bad_tolerances[] = {0, -1, 1, INFINITY, NAN};
    for (size_t i = 0; i < 5; i++) CHECK(m_eigen_symmetric_with_options(&a, bad_tolerances[i], 50, &values, &vectors) == M_ARGUMENT);
    CHECK(m_eigen_symmetric_with_options(&a, 1e-12, 0, &values, &vectors) == M_ARGUMENT);
    CHECK(m_eigen_symmetric_with_options(&a, 1e-12, 10001, &values, &vectors) == M_ARGUMENT);
    CHECK(m_eigen_symmetric(&a, &values, &values) == M_ARGUMENT);
    CHECK(m_eigen_symmetric(&a, &a, &vectors) == M_ARGUMENT);
    CHECK(m_eigen_symmetric(&a, NULL, &vectors) == M_ARGUMENT);
    a.values[0] = NAN; CHECK(m_eigen_symmetric(&a, &values, &vectors) == M_NONFINITE);
    m_free(&a);
    const double asymmetric[] = {1, 2, 2 + 1e-15, 1};
    OK(m_from_array(2, 2, asymmetric, &a)); CHECK(m_eigen_symmetric(&a, &values, &vectors) == M_NOT_SYMMETRIC); m_free(&a);
    OK(m_create(2, 3, &a)); CHECK(m_eigen_symmetric(&a, &values, &vectors) == M_SHAPE); m_free(&a);
    const double overflow[] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};
    OK(m_from_array(2, 2, overflow, &a)); CHECK(m_eigen_symmetric(&a, &values, &vectors) == M_NONFINITE); m_free(&a);
    CHECK(values.values == NULL && vectors.values == NULL && values.cols == 0 && vectors.cols == 0);
    return 0;
}


static int rotation_invariants(const matrix *rotation) {
    const size_t n = rotation->rows;
    CHECK(n == rotation->cols);
    for (size_t row = 0; row < n; row++) for (size_t col = 0; col < n; col++) {
        double dot = 0;
        for (size_t k = 0; k < n; k++) dot += rotation->values[k*n+row] * rotation->values[k*n+col];
        CHECK(fabs(dot - (row == col ? 1.0 : 0.0)) < 2e-14);
    }
    double determinant;
    OK(m_determinant(rotation, &determinant)); CHECK(fabs(determinant - 1) < 2e-14);
    return 0;
}
static int vector_rotation_tests(void) {
    matrix a = {0}, b = {0}, result = {0}, reverse = {0};
    const double av[] = {2,3,4}, bv[] = {5,6,7}, expected[] = {-3,6,-3};
    for (size_t left_rows = 1; left_rows <= 3; left_rows += 2)
    for (size_t right_rows = 1; right_rows <= 3; right_rows += 2) {
        OK(m_from_array(left_rows, 3/left_rows, av, &a));
        OK(m_from_array(right_rows, 3/right_rows, bv, &b));
        OK(m_cross(&a, &b, &result)); OK(m_cross(&b, &a, &reverse));
        CHECK(result.rows == left_rows && result.cols == 3/left_rows);
        double dot_a = 0, dot_b = 0;
        for (size_t i = 0; i < 3; i++) {
            CHECK(result.values[i] == expected[i] && reverse.values[i] == -expected[i]);
            dot_a += av[i]*result.values[i]; dot_b += bv[i]*result.values[i];
        }
        CHECK(dot_a == 0 && dot_b == 0);
        CHECK(memcmp(a.values, av, sizeof(av)) == 0 && memcmp(b.values, bv, sizeof(bv)) == 0);
        m_free(&result); m_free(&reverse);
        OK(m_cross(&a, &a, &result)); CHECK(result.values[0] == 0 && result.values[1] == 0 && result.values[2] == 0);
        m_free(&result); m_free(&a); m_free(&b);
    }
    const double ex[] = {1,0,0}, ey[] = {0,1,0};
    OK(m_from_array(3,1,ex,&a)); OK(m_from_array(3,1,ey,&b)); OK(m_cross(&a,&b,&result));
    CHECK(result.values[0] == 0 && result.values[1] == 0 && result.values[2] == 1);
    m_free(&result);
    CHECK(m_cross(&a,&b,&a) == M_ARGUMENT);
    a.values[0] = INFINITY; CHECK(m_cross(&a,&b,&result) == M_NONFINITE);
    a.values[0] = DBL_MAX; b.values[1] = 2; CHECK(m_cross(&a,&b,&result) == M_NONFINITE);
    CHECK(result.values == NULL && result.rows == 0);
    m_free(&a); m_free(&b);
    OK(m_create(3,3,&a)); OK(m_create(3,1,&b)); CHECK(m_cross(&a,&b,&result) == M_SHAPE);
    m_free(&a); m_free(&b);
    const double quarter = acos(-1.0)/2;
    OK(m_rotation_2d(quarter,&result)); CHECK(fabs(result.values[2]-1) < 1e-15 && fabs(result.values[0]) < 1e-15);
    CHECK(rotation_invariants(&result) == 0); m_free(&result);
    matrix_status (*factories[])(double, matrix *) = {m_rotation_x,m_rotation_y,m_rotation_z};
    const double axes[][3] = {{1,0,0},{0,1,0},{0,0,1}};
    const double angles[] = {0, quarter, -acos(-1.0), 4*acos(-1.0), 1e-8, 0.5, DBL_MAX};
    for (size_t axis = 0; axis < 3; axis++) for (size_t angle = 0; angle < 7; angle++) {
        OK(factories[axis](angles[angle],&result)); CHECK(rotation_invariants(&result) == 0);
        OK(m_from_array(3,1,axes[axis],&a)); OK(m_rotation_axis_angle(&a,angles[angle],&reverse));
        for (size_t i = 0; i < 9; i++) CHECK(fabs(result.values[i]-reverse.values[i]) < 2e-15);
        m_free(&a); m_free(&result); m_free(&reverse);
    }
    OK(m_rotation_x(quarter,&result)); CHECK(fabs(result.values[7]-1) < 1e-15); m_free(&result);
    OK(m_rotation_y(quarter,&result)); CHECK(fabs(result.values[6]+1) < 1e-15); m_free(&result);
    OK(m_rotation_z(quarter,&result)); CHECK(fabs(result.values[3]-1) < 1e-15); m_free(&result);
    const double scales[] = {1,1e-300,1e300,DBL_TRUE_MIN};
    const double normal[] = {1,1,1}; OK(m_from_array(1,3,normal,&a)); OK(m_rotation_axis_angle(&a,.5,&reverse)); m_free(&a);
    for (size_t i = 0; i < 4; i++) {
        const double axis[] = {scales[i],scales[i],scales[i]};
        OK(m_from_array(3,1,axis,&a)); OK(m_rotation_axis_angle(&a,.5,&result)); CHECK(rotation_invariants(&result) == 0);
        for (size_t j = 0; j < 9; j++) CHECK(fabs(result.values[j]-reverse.values[j]) < 2e-15);
        CHECK(memcmp(a.values, axis, sizeof(axis)) == 0);
        for (size_t row = 0; row < 3; row++) CHECK(fabs(result.values[row*3]+result.values[row*3+1]+result.values[row*3+2]-1) < 1e-14);
        m_free(&result); m_free(&a);
    }
    m_free(&reverse);
    const double xy[] = {1,1,0}; OK(m_from_array(3,1,xy,&a)); OK(m_rotation_axis_angle(&a,1e-8,&result));
    CHECK(fabs(result.values[1]/2.5e-17 - 1) < 1e-14); m_free(&result);
    CHECK(m_rotation_axis_angle(&a,NAN,&result) == M_NONFINITE);
    CHECK(m_rotation_axis_angle(&a,.5,&a) == M_ARGUMENT);
    a.values[0] = INFINITY; CHECK(m_rotation_axis_angle(&a,.5,&result) == M_NONFINITE); m_free(&a);
    OK(m_create(3,1,&a)); CHECK(m_rotation_axis_angle(&a,.5,&result) == M_ARGUMENT); m_free(&a);
    OK(m_create(2,1,&a)); CHECK(m_rotation_axis_angle(&a,.5,&result) == M_SHAPE); m_free(&a);
    CHECK(m_rotation_2d(INFINITY,&result) == M_NONFINITE);
    for (size_t i = 0; i < 3; i++) CHECK(factories[i](NAN,&result) == M_NONFINITE);
    CHECK(m_rotation_2d(.5,NULL) == M_ARGUMENT);
    OK(m_identity(2,&result)); CHECK(m_rotation_2d(.5,&result) == M_ARGUMENT && result.values[0] == 1); m_free(&result);
    return 0;
}


static int general_invariants(const matrix *a, const matrix *dr, const matrix *di,
                              const matrix *vr, const matrix *vi) {
    const size_t n = a->rows;
    CHECK(dr->rows == n && dr->cols == 1 && di->rows == n && di->cols == 1);
    CHECK(vr->rows == n && vr->cols == n && vi->rows == n && vi->cols == n);
    double scale = 1;
    for (size_t i = 0; i < n * n; i++) scale = fmax(scale, fabs(a->values[i]));
    for (size_t col = 0; col < n; col++) {
        CHECK(isfinite(dr->values[col]) && isfinite(di->values[col]));
        if (col) CHECK(dr->values[col-1] < dr->values[col] ||
                      (dr->values[col-1] == dr->values[col] && di->values[col-1] <= di->values[col]));
        double norm = 0;
        for (size_t row = 0; row < n; row++) {
            const size_t index = row * n + col;
            double ar = 0, ai = 0;
            for (size_t k = 0; k < n; k++) {
                ar += (a->values[row*n+k]/scale) * vr->values[k*n+col];
                ai += (a->values[row*n+k]/scale) * vi->values[k*n+col];
            }
            const double br = vr->values[index]*(dr->values[col]/scale) - vi->values[index]*(di->values[col]/scale);
            const double bi = vi->values[index]*(dr->values[col]/scale) + vr->values[index]*(di->values[col]/scale);
            CHECK(hypot(ar-br, ai-bi) < 5e-11 * (double)n);
            norm = hypot(norm, hypot(vr->values[index], vi->values[index]));
        }
        CHECK(fabs(norm - 1) < 1e-12 * (double)n);
    }
    return 0;
}

static int general_tests(void) {
    matrix a = {0}, dr = {0}, di = {0}, vr = {0}, vi = {0};
    const double scales[] = {1, 1e300, 1e-300, 0x1p-1074};
    for (size_t k = 0; k < 4; k++) {
        const double input[] = {0, -scales[k], scales[k], 0};
        OK(m_from_array(2, 2, input, &a)); OK(m_eigen_general(&a, &dr, &di, &vr, &vi));
        CHECK(general_invariants(&a,&dr,&di,&vr,&vi) == 0);
        CHECK(dr.values[0] == 0 && dr.values[1] == 0);
        CHECK(fabs(di.values[0]/scales[k]+1)<1e-14 && fabs(di.values[1]/scales[k]-1)<1e-14);
        CHECK(a.values[1] == -scales[k]);
        m_free(&a); m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi);
    }
    for (size_t n = 0; n < 10; n++) {
        OK(m_identity(n,&a));
        for (size_t i = 0; i + 1 < n; i++) a.values[i*n+i+1]=1;
        OK(m_eigen_general(&a,&dr,&di,&vr,&vi)); CHECK(general_invariants(&a,&dr,&di,&vr,&vi)==0);
        for (size_t i = 0; i < n; i++) CHECK(dr.values[i]==1 && di.values[i]==0);
        m_free(&a); m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi);
    }
    uint32_t state = 41;
    for (size_t trial = 0; trial < 40; trial++) {
        OK(m_create(7,7,&a));
        for (size_t i = 0; i < 49; i++) { state=state*1664525u+1013904223u; a.values[i]=(double)(state%101)-50; }
        OK(m_eigen_general(&a,&dr,&di,&vr,&vi)); CHECK(general_invariants(&a,&dr,&di,&vr,&vi)==0);
        m_free(&a); m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi);
    }
    const double diagonal[] = {DBL_MAX,0,0,0,0x1p-1074,0,0,0,-DBL_MAX};
    OK(m_from_array(3,3,diagonal,&a)); OK(m_eigen_general(&a,&dr,&di,&vr,&vi));
    CHECK(dr.values[0]==-DBL_MAX && dr.values[1]==0x1p-1074 && dr.values[2]==DBL_MAX);
    m_free(&a); m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi);
    const double triangles[][4] = {{1e300,1,0,1e-300}, {1e-300,0,1,1e300}};
    for (size_t i = 0; i < 2; i++) {
        OK(m_from_array(2,2,triangles[i],&a)); OK(m_eigen_general(&a,&dr,&di,&vr,&vi));
        CHECK(dr.values[0]==1e-300 && dr.values[1]==1e300);
        CHECK(general_invariants(&a,&dr,&di,&vr,&vi)==0);
        m_free(&a); m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi);
    }
    const double dense[] = {1,2,3,4,5,6,7,8,10};
    OK(m_from_array(3,3,dense,&a));
    CHECK(m_eigen_general_with_options(&a,1,&dr,&di,&vr,&vi)==M_NO_CONVERGENCE);
    CHECK(!dr.values&&!di.values&&!vr.values&&!vi.values);
    CHECK(m_eigen_general_with_options(&a,0,&dr,&di,&vr,&vi)==M_ARGUMENT);
    CHECK(m_eigen_general_with_options(&a,100001,&dr,&di,&vr,&vi)==M_ARGUMENT);
    CHECK(m_eigen_general(&a,&dr,&dr,&vr,&vi)==M_ARGUMENT);
    CHECK(m_eigen_general(&a,&a,&di,&vr,&vi)==M_ARGUMENT);
    CHECK(m_eigen_general(&a,NULL,&di,&vr,&vi)==M_ARGUMENT);
    OK(m_identity(1,&dr)); CHECK(m_eigen_general(&a,&dr,&di,&vr,&vi)==M_ARGUMENT && dr.values[0]==1); m_free(&dr);
    a.values[0]=NAN; CHECK(m_eigen_general(&a,&dr,&di,&vr,&vi)==M_NONFINITE); m_free(&a);
    OK(m_create(2,3,&a)); CHECK(m_eigen_general(&a,&dr,&di,&vr,&vi)==M_SHAPE); m_free(&a);
    const double overflow[] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
    OK(m_from_array(2,2,overflow,&a)); CHECK(m_eigen_general(&a,&dr,&di,&vr,&vi)==M_NONFINITE); m_free(&a);
    const double mixed[] = {0,1e300,-1e-300,0};
    OK(m_from_array(2,2,mixed,&a)); OK(m_eigen_general(&a,&dr,&di,&vr,&vi));
    CHECK(fabs(di.values[0]+1)<1e-14 && fabs(di.values[1]-1)<1e-14);
    CHECK(general_invariants(&a,&dr,&di,&vr,&vi)==0);
    m_free(&a); m_free(&dr); m_free(&di); m_free(&vr); m_free(&vi);
    return 0;
}

static int solver_tests(void) {
    for (int algorithm=M_LU; algorithm<=M_QR; algorithm++) {
        matrix a={0}, b={0}, x={0}; matrix_factor *factor=NULL;
        const double av[]={4,1,1,3}, bv[]={6,5,7,4}, expected[]={1,1,2,1};
        OK(m_from_array(2,2,av,&a)); OK(m_from_array(2,2,bv,&b));
        OK(m_factorize(&a,(matrix_factor_algorithm)algorithm,&factor));
        CHECK(m_factorize(&a,M_LU,&factor)==M_ARGUMENT);
        a.values[0]=99; m_free(&a);
        for(int run=0;run<3;run++) {
            OK(m_factor_solve(factor,&b,&x));
            for(size_t i=0;i<4;i++)CHECK(fabs(x.values[i]-expected[i])<1e-12);
            CHECK(m_factor_solve(factor,&b,&x)==M_ARGUMENT); m_free(&x);
        }
        CHECK(m_factor_solve(factor,&b,&b)==M_ARGUMENT);
        CHECK(b.values[0]==6);
        double rcond=-1; OK(m_factor_rcond(factor,&rcond)); CHECK(fabs(rcond-.44)<1e-12);
        OK(m_create(1,1,&a)); CHECK(m_factor_solve(factor,&a,&x)==M_SHAPE); CHECK(x.values==NULL);
        b.values[0]=INFINITY; CHECK(m_factor_solve(factor,&b,&x)==M_NONFINITE); CHECK(x.values==NULL);
        m_factor_free(&factor); m_factor_free(&factor); CHECK(factor==NULL);
        CHECK(m_factor_solve(factor,&a,&x)==M_ARGUMENT);
        m_free(&a); m_free(&b);
    }
    matrix singular={0}; matrix_factor *factor=NULL;
    const double values[]={1,2,2,4}; OK(m_from_array(2,2,values,&singular));
    CHECK(m_factorize(&singular,M_LU,&factor)==M_SINGULAR); CHECK(factor==NULL);
    CHECK(m_solve(&singular,&singular,&singular)==M_ARGUMENT); CHECK(singular.values[0]==1);
    m_free(&singular);
    return 0;
}

int main(void) {
    CHECK(solver_tests() == 0);
    CHECK(general_tests() == 0);
    CHECK(determinant_tests() == 0);
    CHECK(eigen_tests() == 0);
    CHECK(vector_rotation_tests() == 0);
    double input[] = {1, 2, 3, 4, 5, 6};
    matrix a = {0}, b = {0}, out = {0}, copy = {0};
    OK(m_from_array(2, 3, input, &a));
    input[0] = 100;
    CHECK(a.values[0] == 1);
    OK(m_copy(&a, &copy));
    OK(m_set(&copy, 0, 0, 99));
    CHECK(a.values[0] == 1 && copy.values[0] == 99);
    m_free(&copy);
    OK(m_transpose(&a, &out));
    CHECK(out.rows == 3 && out.cols == 2);
    const double transposed[] = {1, 4, 2, 5, 3, 6};
    for (size_t i = 0; i < 6; i++) CHECK(out.values[i] == transposed[i]);
    m_free(&out);
    OK(m_row(&a, 1, &out));
    CHECK(out.rows == 1 && out.cols == 3 && out.values[0] == 4 && out.values[2] == 6);
    out.values[0] = -10;
    CHECK(a.values[3] == 4);
    m_free(&out);
    OK(m_column(&a, 1, &out));
    CHECK(out.rows == 2 && out.cols == 1 && out.values[0] == 2 && out.values[1] == 5);
    m_free(&out);
    OK(m_add(&a, &a, &out));
    CHECK(out.values[0] == 2 && out.values[5] == 12 && a.values[0] == 1);
    m_free(&out);
    OK(m_subtract(&a, &a, &out));
    for (size_t i = 0; i < 6; i++) CHECK(out.values[i] == 0);
    m_free(&out);
    OK(m_scale(&a, 0.5, &out));
    CHECK(out.values[0] == 0.5 && out.values[5] == 3);
    m_free(&out);
    const double right[] = {7, 8, 9, 10, 11, 12};
    OK(m_from_array(3, 2, right, &b));
    OK(m_multiply(&a, &b, &out));
    const double product[] = {58, 64, 139, 154};
    for (size_t i = 0; i < 4; i++) CHECK(out.values[i] == product[i]);
    m_free(&out);
    CHECK(m_add(&a, &b, &out) == M_SHAPE);
    CHECK(m_multiply(&a, &a, &out) == M_SHAPE);
    CHECK(m_add(&a, &a, &a) == M_ARGUMENT);
    CHECK(m_row(&a, 2, &out) == M_INDEX);
    CHECK(m_column(&a, 3, &out) == M_INDEX);
    double value = -1;
    CHECK(m_get(&a, 2, 0, &value) == M_INDEX);
    CHECK(m_trace(&a, &value) == M_SHAPE && value == -1);
    CHECK(m_determinant(&a, &value) == M_SHAPE);
    CHECK(m_set(&a, 0, 0, INFINITY) == M_NONFINITE && a.values[0] == 1);
    bool upper, lower;
    OK(m_triangular(&a, &upper, &lower)); CHECK(!upper && !lower);
    m_free(&a); m_free(&b);

    const double square[] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
    OK(m_from_array(3, 3, square, &a));
    OK(m_determinant(&a, &value)); CHECK(fabs(value + 306) < 1e-10);
    OK(m_trace(&a, &value)); CHECK(value == 11);
    CHECK(a.values[0] == 6 && a.values[8] == 7);
    m_free(&a);
    const double swapped[] = {0, 1, 1, 0};
    OK(m_from_array(2, 2, swapped, &a));
    OK(m_determinant(&a, &value)); CHECK(value == -1);
    m_free(&a);
    const double singular[] = {1, 2, 2, 4};
    OK(m_from_array(2, 2, singular, &a));
    OK(m_determinant(&a, &value)); CHECK(value == 0);
    m_free(&a);
    OK(m_identity(3, &a));
    OK(m_set(&a, 0, 1, 3));
    OK(m_triangular(&a, &upper, &lower)); CHECK(upper && !lower);
    OK(m_determinant(&a, &value)); CHECK(value == 1);
    m_free(&a);
    OK(m_create(0, 0, &a));
    OK(m_determinant(&a, &value)); CHECK(value == 1);
    OK(m_trace(&a, &value)); CHECK(value == 0);
    m_free(&a);
    OK(m_create(2, 0, &a)); OK(m_create(0, 3, &b));
    OK(m_multiply(&a, &b, &out));
    CHECK(out.rows == 2 && out.cols == 3);
    for (size_t i = 0; i < 6; i++) CHECK(out.values[i] == 0);
    m_free(&out); m_free(&a); m_free(&b);
    CHECK(m_create(SIZE_MAX, 2, &out) == M_MEMORY);
    const double bad[] = {NAN};
    CHECK(m_from_array(1, 1, bad, &out) == M_NONFINITE);
    CHECK(!out.values && !out.rows && !out.cols);
    const double large[] = {1e308};
    OK(m_from_array(1, 1, large, &a));
    CHECK(m_scale(&a, 2, &out) == M_NONFINITE);
    CHECK(!out.values && a.values[0] == 1e308);
    m_free(&a); m_free(&a);
    puts("native matrix tests passed");
    return 0;
}
