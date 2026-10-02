package matrix

import (
	"errors"
	"math"
)

// Factorization owns a snapshot and can solve multiple right-hand sides repeatedly.
// QR uses column pivoting and requires numerical full column rank.
type Factorization struct {
	rows, cols  int
	algorithm   string
	scale, norm float64
	data, tau   []float64
	permutation []int
}

func solverScale(value, scale float64) (float64, error) {
	result := value / scale
	if !isFinite(result) || value != 0 && result == 0 {
		return 0, errors.New("solver scaling exceeds float64 range")
	}
	return result, nil
}

func factorize(source *Matrix, algorithm string) (*Factorization, error) {
	if source == nil {
		return nil, errors.New("nil matrix")
	}
	m, n := source.rows, source.cols
	if m < n || algorithm != "qr" && m != n {
		return nil, errors.New("factorization requires square input, or rows >= columns for QR")
	}
	f := &Factorization{rows: m, cols: n, algorithm: algorithm, scale: 1, data: source.Values(), tau: make([]float64, n), permutation: make([]int, m)}
	for i := range f.permutation {
		f.permutation[i] = i
	}
	maximum := 0.0
	for _, v := range f.data {
		maximum = math.Max(maximum, math.Abs(v))
	}
	if maximum != 0 {
		_, exponent := math.Frexp(maximum)
		f.scale = math.Ldexp(1, exponent-1)
	}
	for i, v := range f.data {
		var err error
		f.data[i], err = solverScale(v, f.scale)
		if err != nil {
			return nil, err
		}
	}
	a := f.data
	for i := 0; i < m; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += math.Abs(a[i*n+j])
		}
		f.norm = math.Max(f.norm, sum)
	}
	switch algorithm {
	case "lu":
		for k := 0; k < n; k++ {
			pivot := k
			for i := k + 1; i < n; i++ {
				if math.Abs(a[i*n+k]) > math.Abs(a[pivot*n+k]) {
					pivot = i
				}
			}
			if a[pivot*n+k] == 0 {
				return nil, errors.New("singular matrix: zero computed LU pivot")
			}
			for j := 0; j < n; j++ {
				a[k*n+j], a[pivot*n+j] = a[pivot*n+j], a[k*n+j]
			}
			f.permutation[k], f.permutation[pivot] = f.permutation[pivot], f.permutation[k]
			for i := k + 1; i < n; i++ {
				a[i*n+k] /= a[k*n+k]
				for j := k + 1; j < n; j++ {
					a[i*n+j] -= a[i*n+k] * a[k*n+j]
				}
			}
		}
	case "cholesky":
		for i := 0; i < n; i++ {
			for j := 0; j < i; j++ {
				if source.values[i*n+j] != source.values[j*n+i] {
					return nil, errors.New("Cholesky requires exact symmetry")
				}
			}
		}
		for i := 0; i < n; i++ {
			for j := 0; j <= i; j++ {
				value := a[i*n+j]
				for k := 0; k < j; k++ {
					value -= a[i*n+k] * a[j*n+k]
				}
				if i == j {
					if value <= 0 {
						return nil, errors.New("Cholesky requires positive computed pivots")
					}
					a[i*n+j] = math.Sqrt(value)
				} else {
					a[i*n+j] = value / a[j*n+j]
				}
			}
		}
	case "qr":
		largest := 0.0
		for j := 0; j < n; j++ {
			norm := 0.0
			for i := 0; i < m; i++ {
				norm = math.Hypot(norm, a[i*n+j])
			}
			largest = math.Max(largest, norm)
		}
		threshold := 2.220446049250313e-16 * float64(m) * largest
		for k := 0; k < n; k++ {
			pivot, norm := k, 0.0
			for j := k; j < n; j++ {
				candidate := 0.0
				for i := k; i < m; i++ {
					candidate = math.Hypot(candidate, a[i*n+j])
				}
				if candidate > norm {
					pivot, norm = j, candidate
				}
			}
			if norm <= threshold {
				return nil, errors.New("QR input is numerically rank deficient")
			}
			for i := 0; i < m; i++ {
				a[i*n+k], a[i*n+pivot] = a[i*n+pivot], a[i*n+k]
			}
			f.permutation[k], f.permutation[pivot] = f.permutation[pivot], f.permutation[k]
			old := a[k*n+k]
			alpha := -math.Copysign(norm, old)
			divisor := old - alpha
			f.tau[k] = (alpha - old) / alpha
			for i := k + 1; i < m; i++ {
				a[i*n+k] /= divisor
			}
			a[k*n+k] = alpha
			for j := k + 1; j < n; j++ {
				dot := a[k*n+j]
				for i := k + 1; i < m; i++ {
					dot += a[i*n+k] * a[i*n+j]
				}
				dot *= f.tau[k]
				a[k*n+j] -= dot
				for i := k + 1; i < m; i++ {
					a[i*n+j] -= a[i*n+k] * dot
				}
			}
		}
	}
	for _, v := range a {
		if !isFinite(v) {
			return nil, errors.New("solver arithmetic exceeds float64 range")
		}
	}
	return f, nil
}

func (m *Matrix) FactorLU() (*Factorization, error)       { return factorize(m, "lu") }
func (m *Matrix) FactorCholesky() (*Factorization, error) { return factorize(m, "cholesky") }
func (m *Matrix) FactorQR() (*Factorization, error)       { return factorize(m, "qr") }
func (m *Matrix) Solve(rhs *Matrix) (*Matrix, error) {
	f, err := m.FactorLU()
	if err != nil {
		return nil, err
	}
	return f.Solve(rhs)
}
func (m *Matrix) LeastSquares(rhs *Matrix) (*Matrix, error) {
	f, err := m.FactorQR()
	if err != nil {
		return nil, err
	}
	return f.Solve(rhs)
}

func (f *Factorization) solve(rhs *Matrix, rescale bool) (*Matrix, error) {
	if f == nil || f.algorithm == "" {
		return nil, errors.New("uninitialized factorization")
	}
	if rhs == nil || rhs.rows != f.rows {
		return nil, errors.New("right-hand side row count must match factorization")
	}
	m, n, p := f.rows, f.cols, rhs.cols
	work := make([]float64, len(rhs.values))
	a := f.data
	for i := 0; i < m; i++ {
		row := i
		if f.algorithm == "lu" {
			row = f.permutation[i]
		}
		for j := 0; j < p; j++ {
			value := rhs.values[row*p+j]
			if rescale {
				var err error
				value, err = solverScale(value, f.scale)
				if err != nil {
					return nil, err
				}
			}
			work[i*p+j] = value
		}
	}
	if f.algorithm == "qr" {
		for k := 0; k < n; k++ {
			for j := 0; j < p; j++ {
				dot := work[k*p+j]
				for i := k + 1; i < m; i++ {
					dot += a[i*n+k] * work[i*p+j]
				}
				dot *= f.tau[k]
				work[k*p+j] -= dot
				for i := k + 1; i < m; i++ {
					work[i*p+j] -= a[i*n+k] * dot
				}
			}
		}
	} else {
		for i := 0; i < n; i++ {
			for j := 0; j < p; j++ {
				value := work[i*p+j]
				for k := 0; k < i; k++ {
					value -= a[i*n+k] * work[k*p+j]
				}
				if f.algorithm == "cholesky" {
					value /= a[i*n+i]
				}
				work[i*p+j] = value
			}
		}
	}
	for i := n - 1; i >= 0; i-- {
		for j := 0; j < p; j++ {
			value := work[i*p+j]
			for k := i + 1; k < n; k++ {
				coefficient := a[i*n+k]
				if f.algorithm == "cholesky" {
					coefficient = a[k*n+i]
				}
				value -= coefficient * work[k*p+j]
			}
			work[i*p+j] = value / a[i*n+i]
		}
	}
	result, _ := Zeros(n, p)
	for i := 0; i < n; i++ {
		row := i
		if f.algorithm == "qr" {
			row = f.permutation[i]
		}
		copy(result.values[row*p:(row+1)*p], work[i*p:(i+1)*p])
	}
	for _, v := range work {
		if !isFinite(v) {
			return nil, errors.New("solver arithmetic exceeds float64 range")
		}
	}
	return result, nil
}
func (f *Factorization) Solve(rhs *Matrix) (*Matrix, error) { return f.solve(rhs, true) }

// ReciprocalCondition computes infinity-norm rcond using an inverse (O(n^3)).
// It is a floating-point diagnostic, not a certified error bound.
func (f *Factorization) ReciprocalCondition() (float64, error) {
	if f == nil || f.algorithm == "" || f.rows != f.cols {
		return 0, errors.New("condition diagnostic requires a square factorization")
	}
	n := f.cols
	if n == 0 {
		return 1, nil
	}
	identity, err := Identity(n)
	if err != nil {
		return 0, err
	}
	inverse, err := f.solve(identity, false)
	if err != nil {
		return 0, nil
	}
	norm := 0.0
	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += math.Abs(inverse.values[i*n+j])
		}
		norm = math.Max(norm, sum)
	}
	return math.Min(1, (1/f.norm)/norm), nil
}
