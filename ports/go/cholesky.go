package matrix

import (
	"errors"
	"math"
	"sort"
)

// SparseCholeskySymbolic owns a reusable pattern. FillSteps uses -1 for originals.
type SparseCholeskySymbolic struct {
	n                                 int
	sourceRP, sourceCI, rp, ci, steps []int
}
type SparseCholesky struct{ lower *CSRMatrix }

func NewSparseCholeskySymbolic(a *CSRMatrix) (*SparseCholeskySymbolic, error) {
	return analyzeCholesky(a, false)
}
func analyzeCholesky(a *CSRMatrix, incomplete bool) (*SparseCholeskySymbolic, error) {
	if a == nil || a.rows != a.cols {
		return nil, errors.New("Cholesky requires square matrix")
	}
	s := &SparseCholeskySymbolic{n: a.rows, sourceRP: append([]int{}, a.rp...), sourceCI: append([]int{}, a.ci...), rp: []int{0}}
	g := make([]map[int]int, s.n)
	for i := range g {
		g[i] = make(map[int]int)
	}
	for i := 0; i < s.n; i++ {
		for p := a.rp[i]; p < a.rp[i+1]; p++ {
			j := a.ci[p]
			if i != j {
				g[i][j] = -1
				g[j][i] = -1
			}
		}
	}
	for k := 0; !incomplete && k < s.n; k++ {
		ns := []int{}
		for j := range g[k] {
			if j > k {
				ns = append(ns, j)
			}
		}
		sort.Ints(ns)
		for u, i := range ns {
			for _, j := range ns[:u] {
				if _, ok := g[i][j]; !ok {
					g[i][j] = k
					g[j][i] = k
				}
			}
		}
	}
	for i := 0; i < s.n; i++ {
		js := []int{}
		for j := range g[i] {
			if j < i {
				js = append(js, j)
			}
		}
		sort.Ints(js)
		for _, j := range js {
			s.ci = append(s.ci, j)
			s.steps = append(s.steps, g[i][j])
		}
		s.ci = append(s.ci, i)
		s.steps = append(s.steps, -1)
		s.rp = append(s.rp, len(s.ci))
	}
	return s, nil
}
func (s *SparseCholeskySymbolic) Size() int { return s.n }
func (s *SparseCholeskySymbolic) NNZ() int  { return len(s.ci) }
func (s *SparseCholeskySymbolic) FillCount() int {
	n := 0
	for _, k := range s.steps {
		if k >= 0 {
			n++
		}
	}
	return n
}
func (s *SparseCholeskySymbolic) RowOffsets() []int    { return append([]int{}, s.rp...) }
func (s *SparseCholeskySymbolic) ColumnIndices() []int { return append([]int{}, s.ci...) }
func (s *SparseCholeskySymbolic) FillSteps() []int     { return append([]int{}, s.steps...) }
func (s *SparseCholeskySymbolic) Factorize(a *CSRMatrix) (*SparseCholesky, error) {
	if s == nil || a == nil || a.rows != s.n || a.cols != s.n || len(a.rp) != len(s.sourceRP) || len(a.ci) != len(s.sourceCI) {
		return nil, errors.New("Cholesky symbolic pattern mismatch")
	}
	for i, v := range a.rp {
		if v != s.sourceRP[i] {
			return nil, errors.New("Cholesky symbolic pattern mismatch")
		}
	}
	for i, v := range a.ci {
		if v != s.sourceCI[i] {
			return nil, errors.New("Cholesky symbolic pattern mismatch")
		}
	}
	rows := make([]map[int]float64, s.n)
	for i := range rows {
		rows[i] = make(map[int]float64)
		for p := a.rp[i]; p < a.rp[i+1]; p++ {
			rows[i][a.ci[p]] = a.values[p]
		}
	}
	for i, row := range rows {
		for j, v := range row {
			if v != rows[j][i] {
				return nil, errors.New("Cholesky requires symmetric values")
			}
		}
	}
	rp, ci := s.rp, s.ci
	v := make([]float64, len(ci))
	for i := 0; i < s.n; i++ {
		for p := rp[i]; p < rp[i+1]; p++ {
			j := ci[p]
			t := rows[i][j]
			u, w := rp[i], rp[j]
			for u < p && w < rp[j+1]-1 {
				if ci[u] == ci[w] {
					t -= v[u] * v[w]
					u++
					w++
				} else if ci[u] < ci[w] {
					u++
				} else {
					w++
				}
			}
			if !isFinite(t) {
				return nil, errors.New("nonfinite Cholesky factor")
			}
			if i == j {
				if t <= 0 {
					return nil, errors.New("nonpositive Cholesky pivot")
				}
				v[p] = math.Sqrt(t)
			} else {
				v[p] = t / v[rp[j+1]-1]
			}
			if !isFinite(v[p]) {
				return nil, errors.New("nonfinite Cholesky factor")
			}
		}
	}
	l, err := NewCSR(s.n, s.n, rp, ci, v)
	if err != nil {
		return nil, err
	}
	return &SparseCholesky{l}, nil
}
func (f *SparseCholesky) Size() int { return f.lower.rows }
func (f *SparseCholesky) NNZ() int  { return f.lower.NNZ() }
func (f *SparseCholesky) Lower() *CSRMatrix {
	a := f.lower
	l, _ := NewCSR(a.rows, a.cols, a.rp, a.ci, a.values)
	return l
}
func (f *SparseCholesky) Solve(b []float64) ([]float64, error) {
	if f == nil || f.lower == nil || len(b) != f.lower.rows {
		return nil, errors.New("invalid Cholesky right-hand side")
	}
	for _, z := range b {
		if !isFinite(z) {
			return nil, errors.New("invalid Cholesky right-hand side")
		}
	}
	x := append([]float64{}, b...)
	a := f.lower
	rp, ci, v := a.rp, a.ci, a.values
	for i := 0; i < a.rows; i++ {
		for p := rp[i]; p < rp[i+1]-1; p++ {
			x[i] -= v[p] * x[ci[p]]
		}
		x[i] /= v[rp[i+1]-1]
		if !isFinite(x[i]) {
			return nil, errors.New("nonfinite Cholesky solve")
		}
	}
	for i := a.rows - 1; i >= 0; i-- {
		x[i] /= v[rp[i+1]-1]
		if !isFinite(x[i]) {
			return nil, errors.New("nonfinite Cholesky solve")
		}
		for p := rp[i]; p < rp[i+1]-1; p++ {
			x[ci[p]] -= v[p] * x[i]
			if !isFinite(x[ci[p]]) {
				return nil, errors.New("nonfinite Cholesky solve")
			}
		}
	}
	return x, nil
}

// IC0 owns a zero-fill incomplete Cholesky factor without shifts or pivoting.
type IC0 struct{ factor *SparseCholesky }

func NewIC0(a *CSRMatrix) (*IC0, error) {
	s, e := analyzeCholesky(a, true)
	if e != nil {
		return nil, e
	}
	f, e := s.Factorize(a)
	if e != nil {
		return nil, e
	}
	return &IC0{f}, nil
}
func (f *IC0) Size() int                            { return f.factor.Size() }
func (f *IC0) NNZ() int                             { return f.factor.NNZ() }
func (f *IC0) Lower() *CSRMatrix                    { return f.factor.Lower() }
func (f *IC0) Apply(b []float64) ([]float64, error) { return f.factor.Solve(b) }
