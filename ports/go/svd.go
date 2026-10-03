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

// SpectralDiagnostics describes the original and retained singular spectra.
type SpectralDiagnostics struct {
	Rank                        int
	ReciprocalCondition         float64
	RetainedReciprocalCondition float64
}

func spectralRank(s []float64, cutoff float64) int {
	rank := 0
	for rank < len(s) && s[rank] > 0 && (cutoff == 0 || s[rank]/s[0] > cutoff) {
		rank++
	}
	return rank
}
func (a *Matrix) defaultCutoff() float64 { return float64(max(a.rows, a.cols)) * 2.220446049250313e-16 }
func (a *Matrix) SpectralDiagnostics() (*SpectralDiagnostics, error) {
	return a.SpectralDiagnosticsWith(a.defaultCutoff())
}
func (a *Matrix) SpectralDiagnosticsWith(cutoff float64) (*SpectralDiagnostics, error) {
	if !isFinite(cutoff) || cutoff < 0 || cutoff > 1 {
		return nil, errors.New("invalid relative cutoff")
	}
	r, err := a.SVD()
	if err != nil {
		return nil, err
	}
	rank := spectralRank(r.Values, cutoff)
	d := &SpectralDiagnostics{Rank: rank}
	if len(r.Values) > 0 && r.Values[0] > 0 {
		d.ReciprocalCondition = r.Values[len(r.Values)-1] / r.Values[0]
	}
	if rank > 0 {
		d.RetainedReciprocalCondition = r.Values[rank-1] / r.Values[0]
	}
	return d, nil
}
func (a *Matrix) Pseudoinverse() (*Matrix, error) { return a.PseudoinverseWith(a.defaultCutoff()) }
func (a *Matrix) PseudoinverseWith(cutoff float64) (*Matrix, error) {
	return a.applyInverse(nil, cutoff, 0)
}
func (a *Matrix) SolveMinimumNorm(rhs *Matrix) (*Matrix, error) {
	return a.SolveMinimumNormWith(rhs, a.defaultCutoff())
}
func (a *Matrix) SolveMinimumNormWith(rhs *Matrix, cutoff float64) (*Matrix, error) {
	if rhs == nil {
		return nil, errors.New("missing right-hand side")
	}
	return a.applyInverse(rhs, cutoff, 0)
}
func (a *Matrix) applyInverse(rhs *Matrix, cutoff float64, lambda float64) (*Matrix, error) {
	if !isFinite(cutoff) || cutoff < 0 || cutoff > 1 {
		return nil, errors.New("invalid relative cutoff")
	}
	if rhs != nil && rhs.rows != a.rows {
		return nil, errors.New("incompatible right-hand side")
	}
	r, err := a.SVD()
	if err != nil {
		return nil, err
	}
	m, n, k := a.rows, a.cols, len(r.Values)
	rank := spectralRank(r.Values, cutoff)
	cols := m
	if rhs != nil {
		cols = rhs.cols
	}
	out, err := Zeros(n, cols)
	if err != nil {
		return nil, err
	}
	for j := 0; j < cols; j++ {
		scale := 1.
		if rhs != nil {
			scale = 0
			for i := 0; i < m; i++ {
				scale = math.Max(scale, math.Abs(rhs.values[i*cols+j]))
			}
		}
		if rhs != nil && rank > 0 && scale != 0 {
			for i := 0; i < m; i++ {
				if rhs.values[i*cols+j] != 0 && rhs.values[i*cols+j]/scale == 0 {
					return nil, errors.New("right-hand side scaling discards an entry")
				}
			}
		}
		for p := 0; p < rank; p++ {
			projection := 0.
			if rhs == nil {
				projection = r.U.values[j*k+p]
			} else if scale != 0 {
				for i := 0; i < m; i++ {
					projection += r.U.values[i*k+p] * (rhs.values[i*cols+j] / scale)
				}
			}
			coefficient := 0.
			if rhs != nil {
				coefficient = quotientProduct(projection, r.Values[p], scale)
				if lambda > 0 {
					coefficient = ridgeProduct(projection, r.Values[p], scale, lambda)
				}
				if !isFinite(coefficient) {
					return nil, errors.New("SVD solve outside float64 range")
				}
			}
			for i := 0; i < n; i++ {
				term := r.Vt.values[p*n+i] * coefficient
				if rhs == nil {
					term = quotientProduct(projection, r.Values[p], r.Vt.values[p*n+i])
				}
				out.values[i*cols+j] += term
				if !isFinite(out.values[i*cols+j]) {
					return nil, errors.New("SVD inverse outside float64 range")
				}
			}
		}
	}
	return out, nil
}

// SolveRidge minimizes ||A X-B||² + lambda ||X||². Zero uses the default SVD cutoff.
func (a *Matrix) SolveRidge(rhs *Matrix, lambda float64) (*Matrix, error) {
	if rhs == nil || !isFinite(lambda) || lambda < 0 {
		return nil, errors.New("expected RHS and finite nonnegative lambda")
	}
	cutoff := 0.
	if lambda == 0 {
		cutoff = a.defaultCutoff()
	}
	return a.applyInverse(rhs, cutoff, lambda)
}
func ridgeProduct(p, s, b, lambda float64) float64 {
	if p == 0 || s == 0 || b == 0 {
		return 0
	}
	root := math.Sqrt(lambda)
	d := math.Max(s, root)
	q, t := s/d, root/d
	pf, pe := math.Frexp(p)
	sf, se := math.Frexp(s)
	bf, be := math.Frexp(b)
	df, de := math.Frexp(d)
	return math.Ldexp(pf*sf*bf/(df*df*(q*q+t*t)), pe+se+be-2*de)
}
