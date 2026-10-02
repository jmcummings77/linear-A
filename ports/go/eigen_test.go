package matrix

import (
	"math"
	"reflect"
	"testing"
)

func verifyEigen(t *testing.T, a *Matrix, result *SymmetricEigenDecomposition) {
	t.Helper()
	n, scale := a.rows, 0.0
	for _, value := range a.values {
		scale = math.Max(scale, math.Abs(value))
	}
	if scale == 0 {
		scale = 1
	}
	q := result.Eigenvectors.values
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			aq, dot := 0.0, 0.0
			for k := 0; k < n; k++ {
				aq += (a.values[i*n+k] / scale) * q[k*n+j]
				dot += q[k*n+i] * q[k*n+j]
			}
			wanted := 0.0
			if i == j {
				wanted = 1
			}
			if math.Abs(aq-q[i*n+j]*(result.Eigenvalues[j]/scale)) > 1e-9 {
				t.Fatal("eigenpair residual")
			}
			if math.Abs(dot-wanted) > 1e-10 {
				t.Fatal("orthogonality residual")
			}
		}
	}
	for i := 1; i < n; i++ {
		if result.Eigenvalues[i] < result.Eigenvalues[i-1] {
			t.Fatal("unsorted spectrum")
		}
	}
}

func TestEigenAnalyticAndDenseResiduals(t *testing.T) {
	for n := 0; n <= 12; n++ {
		a, _ := Zeros(n, n)
		for i := 0; i < n; i++ {
			for j := i; j < n; j++ {
				value := float64((i*17+j*13)%19)/8 - 1
				a.values[i*n+j], a.values[j*n+i] = value, value
			}
		}
		original := a.Values()
		result, err := a.EigenSymmetric()
		if err != nil {
			t.Fatal(err)
		}
		verifyEigen(t, a, result)
		if !reflect.DeepEqual(a.values, original) {
			t.Fatal("mutated input")
		}
	}
	for _, scale := range []float64{1, 1e300, 1e-300, 1e-310} {
		a, _ := New(2, 2, []float64{2 * scale, scale, scale, 2 * scale})
		result, err := a.EigenSymmetric()
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(result.Eigenvalues[0]/scale-1) > 1e-12 || math.Abs(result.Eigenvalues[1]/scale-3) > 1e-12 {
			t.Fatal("analytic spectrum")
		}
		verifyEigen(t, a, result)
	}
}

func TestEigenDiagonalRepeatedAndTightTolerance(t *testing.T) {
	a, _ := New(3, 3, []float64{1e300, 0, 0, 0, 1e-300, 0, 0, 0, -2})
	result, err := a.EigenSymmetric()
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(result.Eigenvalues, []float64{-2, 1e-300, 1e300}) {
		t.Fatal("mixed scale diagonal lost")
	}
	identity, _ := Identity(3)
	repeated, _ := identity.EigenSymmetric()
	verifyEigen(t, identity, repeated)
	tiny, _ := New(2, 2, []float64{1, 1e-200, 1e-200, 2})
	result, err = tiny.EigenSymmetricWith(1e-250, 1)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(result.Eigenvectors.values[1]/1e-200-1) > 1e-14 {
		t.Fatal("offdiagonal norm underflow")
	}
}

func TestEigenValidationAndNonconvergence(t *testing.T) {
	rectangular, _ := Zeros(2, 3)
	if _, err := rectangular.EigenSymmetric(); err == nil {
		t.Fatal("accepted nonsquare")
	}
	asymmetric, _ := New(2, 2, []float64{1, 1, 0, 1})
	if _, err := asymmetric.EigenSymmetric(); err == nil {
		t.Fatal("accepted asymmetric")
	}
	a, _ := New(3, 3, []float64{1, 2, 3, 2, 5, 6, 3, 6, 9})
	for _, tolerance := range []float64{0, -1, 1, math.NaN(), math.Inf(1)} {
		if _, err := a.EigenSymmetricWith(tolerance, 50); err == nil {
			t.Fatal("accepted invalid tolerance")
		}
	}
	for _, limit := range []int{0, -1, 1} {
		if _, err := a.EigenSymmetricWith(1e-12, limit); err == nil {
			t.Fatal("accepted invalid or insufficient sweep limit")
		}
	}
	huge, _ := New(2, 2, []float64{math.MaxFloat64, math.MaxFloat64, math.MaxFloat64, math.MaxFloat64})
	if _, err := huge.EigenSymmetric(); err == nil {
		t.Fatal("accepted infinite eigenvalue")
	}
	var nilMatrix *Matrix
	if _, err := nilMatrix.EigenSymmetric(); err == nil {
		t.Fatal("accepted nil matrix")
	}
}
