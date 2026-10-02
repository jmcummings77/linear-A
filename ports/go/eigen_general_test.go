package matrix

import (
	"math"
	"testing"
)

func verifyGeneral(t *testing.T, a *Matrix, g *GeneralEigenDecomposition) {
	t.Helper()
	n := a.rows
	scale := 0.0
	for _, v := range a.values {
		scale = math.Max(scale, math.Abs(v))
	}
	if scale == 0 {
		scale = 1
	}
	if len(g.EigenvaluesReal) != n || len(g.EigenvaluesImag) != n || g.EigenvectorsReal.rows != n || g.EigenvectorsImag.cols != n {
		t.Fatal("shape")
	}
	for j := 0; j < n; j++ {
		if j > 0 && (g.EigenvaluesReal[j] < g.EigenvaluesReal[j-1] || g.EigenvaluesReal[j] == g.EigenvaluesReal[j-1] && g.EigenvaluesImag[j] < g.EigenvaluesImag[j-1]) {
			t.Fatal("ordering")
		}
		norm := 0.0
		r, im := g.EigenvaluesReal[j]/scale, g.EigenvaluesImag[j]/scale
		for i := 0; i < n; i++ {
			vr, vi := g.EigenvectorsReal.values[i*n+j], g.EigenvectorsImag.values[i*n+j]
			norm += vr*vr + vi*vi
			ar, ai := 0.0, 0.0
			for k := 0; k < n; k++ {
				ar += (a.values[i*n+k] / scale) * g.EigenvectorsReal.values[k*n+j]
				ai += (a.values[i*n+k] / scale) * g.EigenvectorsImag.values[k*n+j]
			}
			if !isFinite(ar) || !isFinite(ai) || math.Hypot(ar-r*vr+im*vi, ai-r*vi-im*vr) > 1e-8*float64(max(1, n)) {
				t.Fatalf("residual column %d row %d", j, i)
			}
		}
		if math.Abs(norm-1) > 1e-10 {
			t.Fatalf("vector norm %g", norm)
		}
	}
}
func TestGeneralComplexAndDegenerate(t *testing.T) {
	cases := [][]float64{{}, {7}, {0, -1, 1, 0}, {2, 1, 0, 2}, {0, 1, 0, 0}, {0, 0, 0, 0}, {1e300, 0, 0, 1e-300}, {0, 1e300, -1e-300, 0}}
	for _, data := range cases {
		n := int(math.Sqrt(float64(len(data))))
		a, _ := New(n, n, data)
		before := a.Values()
		g, err := a.EigenGeneral()
		if err != nil {
			t.Fatal(err)
		}
		verifyGeneral(t, a, g)
		for i, v := range before {
			if a.values[i] != v {
				t.Fatal("input mutation")
			}
		}
		if len(data) == 4 && data[1] == 1e300 {
			if math.Abs(g.EigenvaluesImag[0]+1) > 1e-12 || math.Abs(g.EigenvaluesImag[1]-1) > 1e-12 {
				t.Fatal("mixed scale complex pair")
			}
		}
	}
	for _, scale := range []float64{1e300, 1e-300, 1e-310} {
		a, _ := New(2, 2, []float64{0, -scale, scale, 0})
		g, e := a.EigenGeneral()
		if e != nil {
			t.Fatal(e)
		}
		verifyGeneral(t, a, g)
		if math.Abs(g.EigenvaluesImag[1]/scale-1) > 1e-12 {
			t.Fatal("scaled imaginary value")
		}
	}
}
func TestGeneralDenseResidualsAndLimits(t *testing.T) {
	for n := 1; n <= 12; n++ {
		a, _ := Zeros(n, n)
		for i := range a.values {
			a.values[i] = float64((i*17+n*11)%31-15) / 8
		}
		g, e := a.EigenGeneral()
		if e != nil {
			t.Fatalf("n=%d: %v", n, e)
		}
		verifyGeneral(t, a, g)
	}
	a, _ := New(4, 4, []float64{1, 2, 3, 4, 5, 6, 7, 8, 2, 5, 1, 3, 7, 2, 8, 1})
	if _, e := a.EigenGeneralWith(1); e == nil {
		t.Fatal("expected convergence failure")
	}
	for _, limit := range []int{0, -1, 100001} {
		if _, e := a.EigenGeneralWith(limit); e == nil {
			t.Fatal("invalid limit")
		}
	}
	bad, _ := Zeros(2, 3)
	if _, e := bad.EigenGeneral(); e == nil {
		t.Fatal("nonsquare")
	}
	a.values[0] = math.NaN()
	if _, e := a.EigenGeneral(); e == nil {
		t.Fatal("nonfinite")
	}
	var absent *Matrix
	if _, e := absent.EigenGeneral(); e == nil {
		t.Fatal("nil")
	}
}

func TestGeneralTriangularExtremeDiagonals(t *testing.T) {
	for _, values := range [][]float64{{1e300, 1, 0, 1e-300}, {1e300, 0, 1, 1e-300}} {
		a, _ := New(2, 2, values)
		g, e := a.EigenGeneral()
		if e != nil {
			t.Fatal(e)
		}
		if g.EigenvaluesReal[0] != 1e-300 || g.EigenvaluesReal[1] != 1e300 {
			t.Fatal("triangular diagonal was rounded by work scaling")
		}
		if g.EigenvaluesImag[0] != 0 || g.EigenvaluesImag[1] != 0 {
			t.Fatal("triangular spectrum must be real")
		}
		verifyGeneral(t, a, g)
	}
	a, _ := New(2, 2, []float64{math.MaxFloat64, math.MaxFloat64, math.MaxFloat64, math.MaxFloat64})
	if _, e := a.EigenGeneral(); e == nil {
		t.Fatal("nonfinite eigenvalue must fail")
	}
}
