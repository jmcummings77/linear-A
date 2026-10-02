package matrix

import (
	"math"
	"reflect"
	"testing"
)

func geometryAxis(v ...float64) *Matrix { a, _ := New(3, 1, v); return a }
func geometryNear(t *testing.T, a, b *Matrix) {
	t.Helper()
	if a.rows != b.rows || a.cols != b.cols {
		t.Fatal("shape mismatch")
	}
	for i, x := range a.values {
		if math.Abs(x-b.values[i]) > 3e-14 {
			t.Fatalf("%g != %g", x, b.values[i])
		}
	}
}
func TestCrossOrientationsAndValidation(t *testing.T) {
	for _, ar := range []bool{false, true} {
		for _, br := range []bool{false, true} {
			a, b := geometryAxis(1, 2, 3), geometryAxis(4, 5, 6)
			if ar {
				a = a.Transpose()
			}
			if br {
				b = b.Transpose()
			}
			r, e := a.Cross(b)
			if e != nil {
				t.Fatal(e)
			}
			if r.rows != a.rows || r.cols != a.cols || !reflect.DeepEqual(r.values, []float64{-3, 6, -3}) {
				t.Fatal("cross shape/components")
			}
			if !reflect.DeepEqual(a.values, []float64{1, 2, 3}) {
				t.Fatal("mutated input")
			}
		}
	}
	r, _ := geometryAxis(1, 0, 0).Cross(geometryAxis(0, 1, 0))
	if !reflect.DeepEqual(r.values, []float64{0, 0, 1}) {
		t.Fatal("handedness")
	}
	if _, e := geometryAxis(math.MaxFloat64, 0, 0).Cross(geometryAxis(0, 2, 0)); e == nil {
		t.Fatal("overflow")
	}
	if _, e := geometryAxis(1, 2, 3).Cross(nil); e == nil {
		t.Fatal("nil")
	}
	square, _ := Identity(3)
	if _, e := square.Cross(geometryAxis(1, 2, 3)); e == nil {
		t.Fatal("shape")
	}
}
func TestRotationsHandednessAndOrthogonality(t *testing.T) {
	cases := []struct {
		create        func(float64) (*Matrix, error)
		input, output *Matrix
	}{{RotationX, geometryAxis(0, 1, 0), geometryAxis(0, 0, 1)}, {RotationY, geometryAxis(0, 0, 1), geometryAxis(1, 0, 0)}, {RotationZ, geometryAxis(1, 0, 0), geometryAxis(0, 1, 0)}}
	for _, c := range cases {
		r, _ := c.create(math.Pi / 2)
		actual, _ := r.Multiply(c.input)
		geometryNear(t, actual, c.output)
	}
	for _, angle := range []float64{0, 1e-12, -0.5, math.Pi, math.MaxFloat64} {
		for _, create := range []func(float64) (*Matrix, error){Rotation2D, RotationX, RotationY, RotationZ} {
			r, e := create(angle)
			if e != nil {
				t.Fatal(e)
			}
			actual, _ := r.Transpose().Multiply(r)
			identity, _ := Identity(r.rows)
			geometryNear(t, actual, identity)
			inverse, _ := create(-angle)
			geometryNear(t, inverse, r.Transpose())
			det, _ := r.Determinant()
			if math.Abs(det-1) > 3e-14 {
				t.Fatal("determinant")
			}
		}
	}
}
func TestAxisAngleScalingAndErrors(t *testing.T) {
	v := geometryAxis(1, 2, 3)
	expected, _ := RotationAxisAngle(v, 0.5)
	for _, scale := range []float64{1, 1e300, 1e-300, math.SmallestNonzeroFloat64} {
		a := geometryAxis(scale, 2*scale, 3*scale)
		for _, axis := range []*Matrix{a, a.Transpose()} {
			r, e := RotationAxisAngle(axis, 0.5)
			if e != nil {
				t.Fatal(e)
			}
			geometryNear(t, r, expected)
			fixed, _ := r.Multiply(v)
			geometryNear(t, fixed, v)
			inverse, _ := RotationAxisAngle(axis, -0.5)
			product, _ := inverse.Multiply(r)
			identity, _ := Identity(3)
			geometryNear(t, product, identity)
			zero, _ := RotationAxisAngle(axis, 0)
			geometryNear(t, zero, identity)
		}
	}
	for _, axis := range []*Matrix{nil, geometryAxis(0, 0, 0), {rows: 2, cols: 1, values: []float64{1, 2}}} {
		if _, e := RotationAxisAngle(axis, 0); e == nil {
			t.Fatal("invalid axis")
		}
	}
	for _, angle := range []float64{math.NaN(), math.Inf(1), math.Inf(-1)} {
		if _, e := Rotation2D(angle); e == nil {
			t.Fatal("invalid angle")
		}
		if _, e := RotationAxisAngle(v, angle); e == nil {
			t.Fatal("invalid axis angle")
		}
	}
}
