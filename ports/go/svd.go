package matrix

import (
	"errors"
	"math"
	"sort"
)

// SingularValueDecomposition holds economy factors A=U diag(Values) Vt.
// Factors own their storage. Values are nonnegative and descending.
type SingularValueDecomposition struct {
	U      *Matrix
	Values []float64
	Vt     *Matrix
}

func (a *Matrix) SVD() (*SingularValueDecomposition, error) { return a.SVDWith(1e-12, 100) }

// SVDWith uses cyclic one-sided Jacobi rotations; see ports/SVD.md.
func (a *Matrix) SVDWith(tolerance float64, maxSweeps int) (*SingularValueDecomposition, error) {
	fail := func(message string) (*SingularValueDecomposition, error) { return nil, errors.New(message) }
	if !isFinite(tolerance) || tolerance <= 0 || tolerance >= 1 || maxSweeps < 1 || maxSweeps > 10000 {
		return fail("invalid SVD options")
	}
	m, n := a.rows, a.cols
	if m < n {
		at := a.Transpose()
		r, err := at.SVDWith(tolerance, maxSweeps)
		if err != nil {
			return nil, err
		}
		u := r.Vt.Transpose()
		vt := r.U.Transpose()
		return &SingularValueDecomposition{u, r.Values, vt}, nil
	}
	scale := 0.
	for _, x := range a.values {
		if !isFinite(x) {
			return fail("SVD requires finite input")
		}
		scale = math.Max(scale, math.Abs(x))
	}
	if scale == 0 {
		scale = 1
	}
	b := make([]float64, m*n)
	v := make([]float64, n*n)
	for i, x := range a.values {
		b[i] = x / scale
		if x != 0 && b[i] == 0 {
			return fail("SVD scaling discards an entry")
		}
	}
	for j := 0; j < n; j++ {
		v[j*n+j] = 1
	}
	norm := func(j int) float64 {
		s := 0.
		for i := 0; i < m; i++ {
			s = math.Hypot(s, b[i*n+j])
		}
		return s
	}
	converged := false
	for sweep := 0; sweep <= maxSweeps; sweep++ {
		changed := false
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				np, nq := norm(p), norm(q)
				if np == 0 || nq == 0 {
					continue
				}
				corr := 0.
				for i := 0; i < m; i++ {
					corr += (b[i*n+p] / np) * (b[i*n+q] / nq)
				}
				if math.Abs(corr) <= tolerance {
					continue
				}
				changed = true
				if sweep == maxSweeps {
					continue
				}
				pair := math.Max(np, nq)
				ap, aq := np/pair, nq/pair
				delta, g := aq*aq-ap*ap, 2*ap*aq*corr
				t := math.Copysign(1, g)
				if delta != 0 {
					t = g / (delta + math.Copysign(math.Hypot(delta, g), delta))
				}
				if math.Abs(t) < 2.2250738585072014e-308 {
					small := q
					if np < nq {
						small = p
					}
					for i := 0; i < m; i++ {
						b[i*n+small] = 0
					}
					continue
				}
				c := 1 / math.Hypot(1, t)
				s := c * t
				for i := 0; i < m; i++ {
					x, y := b[i*n+p], b[i*n+q]
					b[i*n+p], b[i*n+q] = c*x-s*y, s*x+c*y
				}
				for i := 0; i < n; i++ {
					x, y := v[i*n+p], v[i*n+q]
					v[i*n+p], v[i*n+q] = c*x-s*y, s*x+c*y
				}
			}
		}
		if !changed {
			converged = true
			break
		}
	}
	if !converged {
		return fail("SVD did not converge")
	}
	norms := make([]float64, n)
	order := make([]int, n)
	for j := 0; j < n; j++ {
		norms[j] = norm(j)
		order[j] = j
	}
	sort.SliceStable(order, func(i, j int) bool { return norms[order[i]] > norms[order[j]] })
	u, vt, values := make([]float64, m*n), make([]float64, n*n), make([]float64, n)
	for j, k := range order {
		values[j] = norms[k] * scale
		if !isFinite(values[j]) || (norms[k] != 0 && values[j] == 0) {
			return fail("singular value outside float64 range")
		}
		for i := 0; i < n; i++ {
			vt[j*n+i] = v[i*n+k]
		}
		if norms[k] != 0 {
			for i := 0; i < m; i++ {
				u[i*n+j] = b[i*n+k] / norms[k]
			}
		} else {
			found := false
			for axis := 0; axis < m; axis++ {
				candidate := make([]float64, m)
				candidate[axis] = 1
				for pass := 0; pass < 2; pass++ {
					for col := 0; col < j; col++ {
						dot := 0.
						for i := 0; i < m; i++ {
							dot += candidate[i] * u[i*n+col]
						}
						for i := 0; i < m; i++ {
							candidate[i] -= dot * u[i*n+col]
						}
					}
				}
				length := 0.
				for _, x := range candidate {
					length = math.Hypot(length, x)
				}
				if length > 0.5/math.Sqrt(float64(m)) {
					for i := 0; i < m; i++ {
						u[i*n+j] = candidate[i] / length
					}
					found = true
					break
				}
			}
			if !found {
				return fail("cannot complete SVD null basis")
			}
		}
	}
	um, err := New(m, n, u)
	if err != nil {
		return nil, err
	}
	vm, err := New(n, n, vt)
	if err != nil {
		return nil, err
	}
	return &SingularValueDecomposition{um, values, vm}, nil
}
