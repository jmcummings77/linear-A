package matrix

import (
	"errors"
	"math"
)

// DeterminantAlgorithm selects a floating-point determinant method.
type DeterminantAlgorithm int

const (
	// DeterminantAuto uses an exact-zero triangular shortcut, otherwise LU.
	DeterminantAuto DeterminantAlgorithm = iota
	DeterminantLU
	// DeterminantCholesky requires exactly symmetric, positive-definite input.
	DeterminantCholesky
)

type scaledProduct struct {
	mantissa float64
	exponent int64
	zero     bool
}

func (p *scaledProduct) multiply(value float64) {
	if value == 0 {
		p.zero = true
		return
	}
	fraction, exponent := math.Frexp(value)
	p.mantissa *= fraction * 2
	p.exponent += int64(exponent - 1)
	if math.Abs(p.mantissa) >= 2 {
		p.mantissa *= 0.5
		p.exponent++
	}
}

func (p scaledProduct) result() (float64, error) {
	if p.zero {
		return 0, nil
	}
	if p.exponent > 1023 {
		return 0, errors.New("matrix determinant overflowed")
	}
	if p.exponent < -1075 {
		return math.Copysign(0, p.mantissa), nil
	}
	value := math.Ldexp(p.mantissa, int(p.exponent))
	if !isFinite(value) {
		return 0, errors.New("matrix determinant overflowed")
	}
	return value, nil
}

// quotientProduct retains a representable product when numerator/pivot alone
// would underflow. The ordinary elimination path avoids this extra work.
func quotientProduct(numerator, pivot, value float64) float64 {
	if value == 0 {
		return 0
	}
	a, ae := math.Frexp(numerator)
	b, be := math.Frexp(pivot)
	c, ce := math.Frexp(value)
	return math.Ldexp(float64(a/b*c), ae-be+ce)
}

// Determinant computes det(0x0)=1 and uses Auto without changing the input.
func (m *Matrix) Determinant() (float64, error) {
	return m.DeterminantWith(DeterminantAuto)
}

// DeterminantWith uses O(n^3) elimination and O(n^2) working storage, except
// triangular Auto inputs. Partial pivoting and reversible binary row scaling
// reduce avoidable overflow/underflow; rounding and conditioning still apply.
// Cholesky rejects nonpositive computed pivots, including rounding failures.
func (m *Matrix) DeterminantWith(algorithm DeterminantAlgorithm) (float64, error) {
	if algorithm < DeterminantAuto || algorithm > DeterminantCholesky {
		return 0, errors.New("unknown determinant algorithm")
	}
	if err := m.square(); err != nil {
		return 0, err
	}
	for _, value := range m.values {
		if !isFinite(value) {
			return 0, errors.New("matrix values must be finite")
		}
	}
	if algorithm == DeterminantCholesky {
		return m.choleskyDeterminant()
	}
	if algorithm == DeterminantAuto {
		upper, lower := m.Triangular()
		if upper || lower {
			product := scaledProduct{mantissa: 1}
			for i := 0; i < m.rows; i++ {
				product.multiply(m.values[i*m.cols+i])
			}
			return product.result()
		}
	}
	size, work := m.rows, m.Values()
	product := scaledProduct{mantissa: 1}
	for row := 0; row < size; row++ {
		largest := 0.0
		for col := 0; col < size; col++ {
			largest = math.Max(largest, math.Abs(work[row*size+col]))
		}
		if largest == 0 {
			return 0, nil
		}
		_, exponent := math.Frexp(largest)
		exponent--
		preserves := true
		for col := 0; col < size; col++ {
			value := work[row*size+col]
			scaled := math.Ldexp(value, -exponent)
			if math.Ldexp(scaled, exponent) != value {
				preserves = false
				break
			}
		}
		if preserves {
			for col := 0; col < size; col++ {
				work[row*size+col] = math.Ldexp(work[row*size+col], -exponent)
			}
			product.exponent += int64(exponent)
		}
	}
	for col := 0; col < size; col++ {
		pivotRow := col
		for row := col + 1; row < size; row++ {
			if math.Abs(work[row*size+col]) > math.Abs(work[pivotRow*size+col]) {
				pivotRow = row
			}
		}
		if work[pivotRow*size+col] == 0 {
			return 0, nil
		}
		if pivotRow != col {
			for j := col; j < size; j++ {
				work[col*size+j], work[pivotRow*size+j] = work[pivotRow*size+j], work[col*size+j]
			}
			product.mantissa = -product.mantissa
		}
		pivot := work[col*size+col]
		product.multiply(pivot)
		for row := col + 1; row < size; row++ {
			numerator := work[row*size+col]
			factor := numerator / pivot
			scaledUpdate := numerator != 0 && math.Abs(factor) < 0x1p-1022
			work[row*size+col] = 0
			for j := col + 1; j < size; j++ {
				term := float64(factor * work[col*size+j])
				if scaledUpdate {
					term = quotientProduct(numerator, pivot, work[col*size+j])
				}
				work[row*size+j] -= term
				if !isFinite(work[row*size+j]) {
					return 0, errors.New("matrix elimination overflowed")
				}
			}
		}
	}
	return product.result()
}

func (m *Matrix) choleskyDeterminant() (float64, error) {
	size, work := m.rows, m.Values()
	for row := 0; row < size; row++ {
		for col := 0; col < row; col++ {
			if work[row*size+col] != work[col*size+row] {
				return 0, errors.New("Cholesky requires exact symmetry")
			}
		}
	}
	product := scaledProduct{mantissa: 1}
	for row := 0; row < size; row++ {
		for col := 0; col <= row; col++ {
			sum := work[row*size+col]
			for k := 0; k < col; k++ {
				term := float64(work[row*size+k] * work[col*size+k])
				sum -= term
			}
			if !isFinite(sum) || row == col && sum <= 0 {
				return 0, errors.New("Cholesky requires finite positive pivots")
			}
			if row == col {
				work[row*size+col] = math.Sqrt(sum)
			} else {
				work[row*size+col] = sum / work[col*size+col]
			}
			if !isFinite(work[row*size+col]) {
				return 0, errors.New("Cholesky produced a nonfinite factor")
			}
		}
		product.multiply(work[row*size+row])
		product.multiply(work[row*size+row])
	}
	return product.result()
}
