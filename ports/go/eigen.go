package matrix

import (
	"errors"
	"math"
	"sort"
)

// SymmetricEigenDecomposition contains ascending eigenvalues and corresponding
// orthonormal eigenvector columns. Repeated eigenvalues can have any orthonormal
// basis of their eigenspace. The result owns its storage.
type SymmetricEigenDecomposition struct {
	Eigenvalues  []float64
	Eigenvectors *Matrix
}

// EigenSymmetric solves a finite, exactly symmetric square matrix with cyclic
// Jacobi rotations, relative Frobenius tolerance 1e-12 and at most 50 sweeps.
// It does not modify its input. Tiny eigenvalues in mixed-scale matrices can
// have large relative errors even when the normwise residual is small.
func (m *Matrix) EigenSymmetric() (*SymmetricEigenDecomposition, error) {
	return m.EigenSymmetricWith(1e-12, 50)
}

// EigenSymmetricWith sets the relative tolerance (finite, strictly between zero
// and one) and positive sweep limit. Nonconvergence and nonfinite results fail.
func (m *Matrix) EigenSymmetricWith(tolerance float64, maxSweeps int) (*SymmetricEigenDecomposition, error) {
	if m == nil || m.rows != m.cols {
		return nil, errors.New("eigendecomposition requires a square matrix")
	}
	if !isFinite(tolerance) || tolerance <= 0 || tolerance >= 1 || maxSweeps <= 0 {
		return nil, errors.New("invalid eigensolver tolerance or sweep limit")
	}
	n := m.rows
	scale := 0.0
	diagonal := true
	for row := 0; row < n; row++ {
		for col := 0; col < n; col++ {
			value := m.values[row*n+col]
			if !isFinite(value) || value != m.values[col*n+row] {
				return nil, errors.New("eigendecomposition requires a finite exactly symmetric matrix")
			}
			scale = math.Max(scale, math.Abs(value))
			if row != col && value != 0 {
				diagonal = false
			}
		}
	}
	vectors, err := Identity(n)
	if err != nil {
		return nil, err
	}
	work := m.Values()
	if diagonal {
		return finishEigen(work, vectors, 1)
	}
	exponentBits := (math.Float64bits(scale) >> 52) & 0x7ff
	if exponentBits != 0 {
		scale = math.Float64frombits(exponentBits << 52)
	}
	norm := 0.0
	for i := range work {
		work[i] /= scale
		norm = math.Hypot(norm, work[i])
	}
	threshold := tolerance * norm
	pairThreshold := threshold / (2 * float64(n))
	for sweep := 0; ; sweep++ {
		offNorm := 0.0
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				offNorm = math.Hypot(offNorm, math.Sqrt2*work[p*n+q])
			}
		}
		if offNorm <= threshold {
			return finishEigen(work, vectors, scale)
		}
		if sweep == maxSweeps {
			return nil, errors.New("symmetric eigendecomposition did not converge within the sweep limit")
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				b := work[p*n+q]
				if math.Abs(b) <= pairThreshold {
					continue
				}
				delta := (work[q*n+q] - work[p*n+p]) / 2
				hypotenuse := math.Hypot(delta, b)
				if delta < 0 {
					hypotenuse = -hypotenuse
				}
				t := b / (delta + hypotenuse)
				c := 1 / math.Sqrt(1+t*t)
				s := t * c
				work[p*n+p] -= t * b
				work[q*n+q] += t * b
				work[p*n+q], work[q*n+p] = 0, 0
				for k := 0; k < n; k++ {
					if k != p && k != q {
						x, y := work[k*n+p], work[k*n+q]
						xp, yq := c*x-s*y, s*x+c*y
						work[k*n+p], work[p*n+k] = xp, xp
						work[k*n+q], work[q*n+k] = yq, yq
					}
					vp, vq := vectors.values[k*n+p], vectors.values[k*n+q]
					vectors.values[k*n+p], vectors.values[k*n+q] = c*vp-s*vq, s*vp+c*vq
				}
			}
		}
	}
}

func finishEigen(work []float64, vectors *Matrix, scale float64) (*SymmetricEigenDecomposition, error) {
	n := vectors.rows
	values, order := make([]float64, n), make([]int, n)
	for i := range values {
		values[i] = work[i*n+i] * scale
		if !isFinite(values[i]) {
			return nil, errors.New("eigenvalue is outside the finite float64 range")
		}
		order[i] = i
	}
	sort.SliceStable(order, func(i, j int) bool { return values[order[i]] < values[order[j]] })
	sorted, err := Zeros(n, n)
	if err != nil {
		return nil, err
	}
	sortedValues := make([]float64, n)
	for col, original := range order {
		sortedValues[col] = values[original]
		length, largest := 0.0, 0
		for row := 0; row < n; row++ {
			value := vectors.values[row*n+original]
			length = math.Hypot(length, value)
			if math.Abs(value) > math.Abs(vectors.values[largest*n+original]) {
				largest = row
			}
		}
		if vectors.values[largest*n+original] < 0 {
			length = -length
		}
		for row := 0; row < n; row++ {
			sorted.values[row*n+col] = vectors.values[row*n+original] / length
		}
	}
	return &SymmetricEigenDecomposition{Eigenvalues: sortedValues, Eigenvectors: sorted}, nil
}
