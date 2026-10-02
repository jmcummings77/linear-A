package matrix

import (
	"errors"
	"math"
)

func vector3(a *Matrix) bool {
	return a != nil && ((a.rows == 3 && a.cols == 1) || (a.rows == 1 && a.cols == 3))
}

// Cross returns the right-handed cross product of two 3x1 or 1x3 vectors,
// preserving the left shape without changing either input. Overflow fails.
func (a *Matrix) Cross(b *Matrix) (*Matrix, error) {
	if !vector3(a) || !vector3(b) {
		return nil, errors.New("cross product requires 3x1 or 1x3 vectors")
	}
	result, err := Zeros(a.rows, a.cols)
	if err != nil {
		return nil, err
	}
	for i := 0; i < 3; i++ {
		j, k := (i+1)%3, (i+2)%3
		value := a.values[j]*b.values[k] - a.values[k]*b.values[j]
		if !isFinite(value) {
			return nil, errors.New("cross product produced a nonfinite component")
		}
		result.values[i] = value
	}
	return result, nil
}

// Rotation2D creates a counterclockwise active column-vector rotation in radians.
func Rotation2D(radians float64) (*Matrix, error) {
	if !isFinite(radians) {
		return nil, errors.New("rotation angle must be finite")
	}
	s, c := math.Sincos(radians)
	return New(2, 2, []float64{c, -s, s, c})
}

// RotationX creates a right-handed active rotation about positive x, in radians.
func RotationX(radians float64) (*Matrix, error) { return rodrigues(1, 0, 0, radians) }

// RotationY creates a right-handed active rotation about positive y, in radians.
func RotationY(radians float64) (*Matrix, error) { return rodrigues(0, 1, 0, radians) }

// RotationZ creates a right-handed active rotation about positive z, in radians.
func RotationZ(radians float64) (*Matrix, error) { return rodrigues(0, 0, 1, radians) }

// RotationAxisAngle creates a right-handed active column-vector rotation in
// radians about a finite nonzero 3x1 or 1x3 axis, normalized without modifying it.
func RotationAxisAngle(axis *Matrix, radians float64) (*Matrix, error) {
	if !vector3(axis) {
		return nil, errors.New("rotation axis requires a 3x1 or 1x3 vector")
	}
	scale := 0.0
	for _, value := range axis.values {
		if !isFinite(value) {
			return nil, errors.New("rotation axis must be finite")
		}
		scale = math.Max(scale, math.Abs(value))
	}
	if scale == 0 {
		return nil, errors.New("rotation axis must be nonzero")
	}
	x, y, z := axis.values[0]/scale, axis.values[1]/scale, axis.values[2]/scale
	norm := math.Sqrt(x*x + y*y + z*z)
	return rodrigues(x/norm, y/norm, z/norm, radians)
}
func rodrigues(x, y, z, radians float64) (*Matrix, error) {
	if !isFinite(radians) {
		return nil, errors.New("rotation angle must be finite")
	}
	s, c := math.Sincos(radians)
	t := 1 - c
	if math.Abs(radians) < 1 {
		half := math.Sin(radians / 2)
		t = 2 * half * half
	}
	return New(3, 3, []float64{c + x*x*t, x*y*t - z*s, x*z*t + y*s, y*x*t + z*s, c + y*y*t, y*z*t - x*s, z*x*t - y*s, z*y*t + x*s, c + z*z*t})
}
