package matrix

import (
	"math"
	"reflect"
	"testing"
)

func mustMatrix(t *testing.T, rows, cols int, values []float64) *Matrix {
	t.Helper()
	m, err := New(rows, cols, values)
	if err != nil {
		t.Fatal(err)
	}
	return m
}

func assertValues(t *testing.T, m *Matrix, expected []float64) {
	t.Helper()
	if !reflect.DeepEqual(m.Values(), expected) {
		t.Fatalf("values = %v, want %v", m.Values(), expected)
	}
}

func TestArithmeticAndRectangularProduct(t *testing.T) {
	a := mustMatrix(t, 2, 3, []float64{1, 2, 3, 4, 5, 6})
	b := mustMatrix(t, 3, 2, []float64{7, 8, 9, 10, 11, 12})
	product, err := a.Multiply(b)
	if err != nil {
		t.Fatal(err)
	}
	assertValues(t, product, []float64{58, 64, 139, 154})
	added, err := a.Add(a)
	if err != nil {
		t.Fatal(err)
	}
	assertValues(t, added, []float64{2, 4, 6, 8, 10, 12})
	scaled, err := a.Scale(2)
	if err != nil {
		t.Fatal(err)
	}
	assertValues(t, scaled, added.Values())
	subtracted, err := a.Subtract(a)
	if err != nil {
		t.Fatal(err)
	}
	assertValues(t, subtracted, []float64{0, 0, 0, 0, 0, 0})
	assertValues(t, a.Transpose(), []float64{1, 4, 2, 5, 3, 6})
	assertValues(t, a.Transpose().Transpose(), a.Values())
	identity, err := Identity(3)
	if err != nil {
		t.Fatal(err)
	}
	unchanged, err := a.Multiply(identity)
	if err != nil {
		t.Fatal(err)
	}
	assertValues(t, unchanged, a.Values())
}

func TestCopyAndIndexingOwnStorage(t *testing.T) {
	source := []float64{1, 2, 3, 4}
	a := mustMatrix(t, 2, 2, source)
	source[0] = 99
	copy := a.Copy()
	if err := copy.Set(0, 0, 42); err != nil {
		t.Fatal(err)
	}
	value, err := a.At(0, 0)
	if err != nil || value != 1 {
		t.Fatalf("At = %g, %v", value, err)
	}
	row, err := a.Row(1)
	if err != nil || !reflect.DeepEqual(row, []float64{3, 4}) {
		t.Fatalf("Row = %v, %v", row, err)
	}
	column, err := a.Column(1)
	if err != nil || !reflect.DeepEqual(column, []float64{2, 4}) {
		t.Fatalf("Column = %v, %v", column, err)
	}
	row[0], column[0] = 99, 99
	values := a.Values()
	values[0] = 99
	assertValues(t, a, []float64{1, 2, 3, 4})
	if _, err := a.At(-1, 0); err == nil {
		t.Fatal("negative index accepted")
	}
	if _, err := a.At(2, 0); err == nil {
		t.Fatal("invalid index accepted")
	}
	if err := a.Set(0, 2, 1); err == nil {
		t.Fatal("invalid index accepted")
	}
	if _, err := a.Row(2); err == nil {
		t.Fatal("invalid row accepted")
	}
	if _, err := a.Column(2); err == nil {
		t.Fatal("invalid column accepted")
	}
}

func TestDeterminantPivotingSingularityAndSign(t *testing.T) {
	cases := []struct {
		size   int
		values []float64
		want   float64
	}{
		{2, []float64{0, 2, 3, 4}, -6},
		{3, []float64{6, 1, 1, 4, -2, 5, 2, 8, 7}, -306},
		{2, []float64{1, 2, 2, 4}, 0},
		{1, []float64{-7}, -7},
		{2, []float64{1e-200, 0, 0, 1}, 1e-200},
	}
	for _, tc := range cases {
		a := mustMatrix(t, tc.size, tc.size, tc.values)
		got, err := a.Determinant()
		if err != nil || math.Abs(got-tc.want) > math.Abs(tc.want)*1e-14 {
			t.Fatalf("det(%v) = %g, %v; want %g", tc.values, got, err, tc.want)
		}
		assertValues(t, a, tc.values)
	}
	trace, err := mustMatrix(t, 2, 2, []float64{0, 2, 3, 4}).Trace()
	if err != nil || trace != 4 {
		t.Fatalf("trace = %g, %v", trace, err)
	}
}

func TestEmptyShapesAndTriangular(t *testing.T) {
	var empty Matrix
	trace, err := empty.Trace()
	if err != nil || trace != 0 {
		t.Fatalf("empty trace = %g, %v", trace, err)
	}
	det, err := empty.Determinant()
	if err != nil || det != 1 {
		t.Fatalf("empty determinant = %g, %v", det, err)
	}
	cases := []struct {
		matrix       *Matrix
		upper, lower bool
	}{
		{&empty, true, true},
		{mustMatrix(t, 2, 2, []float64{1, 0, 0, 1}), true, true},
		{mustMatrix(t, 2, 2, []float64{1, 2, 0, 3}), true, false},
		{mustMatrix(t, 2, 2, []float64{1, 0, 2, 3}), false, true},
		{mustMatrix(t, 2, 0, nil), false, false},
	}
	for _, tc := range cases {
		upper, lower := tc.matrix.Triangular()
		if upper != tc.upper || lower != tc.lower {
			t.Fatalf("triangular = %v, %v; want %v, %v", upper, lower, tc.upper, tc.lower)
		}
	}
	a, _ := Zeros(2, 0)
	b, _ := Zeros(0, 3)
	product, err := a.Multiply(b)
	if err != nil {
		t.Fatal(err)
	}
	assertValues(t, product, []float64{0, 0, 0, 0, 0, 0})
	row, err := a.Row(1)
	if err != nil || len(row) != 0 {
		t.Fatalf("empty row = %v, %v", row, err)
	}
	column, err := b.Column(1)
	if err != nil || len(column) != 0 {
		t.Fatalf("empty column = %v, %v", column, err)
	}
}

func TestInvalidShapesValuesAndOverflow(t *testing.T) {
	if _, err := New(2, 2, []float64{1}); err == nil {
		t.Fatal("wrong value count accepted")
	}
	if _, err := Zeros(-1, 2); err == nil {
		t.Fatal("negative dimensions accepted")
	}
	if _, err := Zeros(int(^uint(0)>>1), 2); err == nil {
		t.Fatal("overflow dimensions accepted")
	}
	if _, err := New(1, 1, []float64{math.NaN()}); err == nil {
		t.Fatal("NaN accepted")
	}
	rectangular, _ := Zeros(2, 3)
	square, _ := Identity(2)
	if _, err := rectangular.Add(square); err == nil {
		t.Fatal("mismatched addition accepted")
	}
	if _, err := rectangular.Subtract(square); err == nil {
		t.Fatal("mismatched subtraction accepted")
	}
	if _, err := rectangular.Multiply(square); err == nil {
		t.Fatal("mismatched multiplication accepted")
	}
	if _, err := rectangular.Trace(); err == nil {
		t.Fatal("nonsquare trace accepted")
	}
	if _, err := rectangular.Determinant(); err == nil {
		t.Fatal("nonsquare determinant accepted")
	}
	if _, err := square.Scale(math.Inf(1)); err == nil {
		t.Fatal("infinite scale accepted")
	}
	if err := square.Set(0, 0, math.NaN()); err == nil {
		t.Fatal("NaN mutation accepted")
	}
	huge := mustMatrix(t, 1, 1, []float64{math.MaxFloat64})
	if _, err := huge.Add(huge); err == nil {
		t.Fatal("addition overflow accepted")
	}
	if _, err := huge.Scale(2); err == nil {
		t.Fatal("scale overflow accepted")
	}
	if _, err := huge.Multiply(huge); err == nil {
		t.Fatal("multiplication overflow accepted")
	}
	if _, err := mustMatrix(t, 2, 2, []float64{math.MaxFloat64, 0, 0, math.MaxFloat64}).Trace(); err == nil {
		t.Fatal("trace overflow accepted")
	}
	if _, err := mustMatrix(t, 2, 2, []float64{math.MaxFloat64, 0, 0, 2}).Determinant(); err == nil {
		t.Fatal("determinant overflow accepted")
	}
}
