package matrix

import (
	"errors"
	"math"
	"sort"
)

// CSRMatrix owns canonical CSR arrays. Accessors return copies.
type CSRMatrix struct {
	rows, cols int
	rp, ci     []int
	values     []float64
}
type CGOptions struct {
	RelativeTolerance, AbsoluteTolerance float64
	MaxIterations                        int
	Jacobi, Capture                      bool
}

func DefaultCGOptions() CGOptions { return CGOptions{RelativeTolerance: 1e-10, MaxIterations: 1000} }

type CGResult struct {
	X          []float64
	Converged  bool
	Iterations int
	Reason     string
	Residuals  []float64
	Iterates   [][]float64
}

func NewCSR(rows, cols int, rp, ci []int, values []float64) (*CSRMatrix, error) {
	if rows < 0 || cols < 0 || rows == int(^uint(0)>>1) || len(rp) != rows+1 || len(ci) != len(values) || rp[0] != 0 || rp[rows] != len(values) {
		return nil, errors.New("invalid CSR dimensions or arrays")
	}
	for _, v := range rp {
		if v < 0 || v > len(values) {
			return nil, errors.New("invalid CSR offsets")
		}
	}
	for i := 0; i < rows; i++ {
		if rp[i] > rp[i+1] {
			return nil, errors.New("CSR offsets must be monotone")
		}
		previous := -1
		for p := rp[i]; p < rp[i+1]; p++ {
			if ci[p] <= previous || ci[p] >= cols || !isFinite(values[p]) {
				return nil, errors.New("CSR columns must be sorted unique and finite")
			}
			previous = ci[p]
		}
	}
	return &CSRMatrix{rows, cols, append([]int{}, rp...), append([]int{}, ci...), append([]float64{}, values...)}, nil
}
func CSRFromDense(a *Matrix) (*CSRMatrix, error) {
	rp := []int{0}
	ci := []int{}
	v := []float64{}
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			x := a.values[i*a.cols+j]
			if x != 0 {
				ci = append(ci, j)
				v = append(v, x)
			}
		}
		rp = append(rp, len(v))
	}
	return NewCSR(a.rows, a.cols, rp, ci, v)
}
func (a *CSRMatrix) Rows() int            { return a.rows }
func (a *CSRMatrix) Cols() int            { return a.cols }
func (a *CSRMatrix) NNZ() int             { return len(a.values) }
func (a *CSRMatrix) RowOffsets() []int    { return append([]int{}, a.rp...) }
func (a *CSRMatrix) ColumnIndices() []int { return append([]int{}, a.ci...) }
func (a *CSRMatrix) Values() []float64    { return append([]float64{}, a.values...) }
func (a *CSRMatrix) Matvec(x []float64) ([]float64, error) {
	if len(x) != a.cols {
		return nil, errors.New("invalid vector")
	}
	for _, v := range x {
		if !isFinite(v) {
			return nil, errors.New("nonfinite vector")
		}
	}
	out := make([]float64, a.rows)
	for i := range out {
		for p := a.rp[i]; p < a.rp[i+1]; p++ {
			out[i] += a.values[p] * x[a.ci[p]]
		}
		if !isFinite(out[i]) {
			return nil, errors.New("sparse multiplication outside float64 range")
		}
	}
	return out, nil
}
func (a *CSRMatrix) find(row, col int) int {
	return a.rp[row] + sort.Search(a.rp[row+1]-a.rp[row], func(k int) bool { return a.ci[a.rp[row]+k] >= col })
}
func (a *CSRMatrix) ConjugateGradient(b []float64, o CGOptions) (*CGResult, error) {
	n := a.rows
	if a.cols != n || len(b) != n {
		return nil, errors.New("CG requires square matrix and matching vector")
	}
	for _, x := range b {
		if !isFinite(x) {
			return nil, errors.New("nonfinite RHS")
		}
	}
	if !isFinite(o.RelativeTolerance) || o.RelativeTolerance < 0 || o.RelativeTolerance >= 1 || !isFinite(o.AbsoluteTolerance) || o.AbsoluteTolerance < 0 || o.MaxIterations < 0 || o.MaxIterations > 100000 {
		return nil, errors.New("invalid CG options")
	}
	diag := make([]float64, n)
	for i := 0; i < n; i++ {
		diag[i] = 1
		for p := a.rp[i]; p < a.rp[i+1]; p++ {
			j := a.ci[p]
			q := a.find(j, i)
			other := 0.
			if q < a.rp[j+1] && a.ci[q] == i {
				other = a.values[q]
			}
			if a.values[p] != other {
				return nil, errors.New("CG requires exact symmetry")
			}
		}
		if o.Jacobi {
			q := a.find(i, i)
			if q == a.rp[i+1] || a.ci[q] != i || a.values[q] <= 0 {
				return nil, errors.New("Jacobi requires positive diagonal")
			}
			diag[i] = a.values[q]
		}
	}
	norm := func(v []float64) float64 {
		s := 0.
		for _, x := range v {
			s = math.Hypot(s, x)
		}
		return s
	}
	dot := func(a, b []float64) float64 {
		s := 0.
		for i, x := range a {
			s += x * b[i]
		}
		return s
	}
	x := make([]float64, n)
	r := append([]float64{}, b...)
	history := []float64{norm(r)}
	frames := [][]float64{}
	if o.Capture {
		frames = append(frames, append([]float64{}, x...))
	}
	threshold := math.Max(o.AbsoluteTolerance, o.RelativeTolerance*history[0])
	result := func(reason string) (*CGResult, error) {
		return &CGResult{append([]float64{}, x...), reason == "converged", len(history) - 1, reason, history, frames}, nil
	}
	if !isFinite(history[0]) {
		return result("nonfinite")
	}
	if history[0] <= threshold {
		return result("converged")
	}
	z := make([]float64, n)
	for i := range z {
		z[i] = r[i] / diag[i]
	}
	p := append([]float64{}, z...)
	rho := dot(r, z)
	for step := 0; step < o.MaxIterations; step++ {
		q, err := a.Matvec(p)
		if err != nil {
			return result("nonfinite")
		}
		curvature := dot(p, q)
		if !isFinite(rho) || !isFinite(curvature) {
			return result("nonfinite")
		}
		if rho <= 0 || curvature <= 0 {
			return result("breakdown")
		}
		alpha := rho / curvature
		candidate := make([]float64, n)
		for i := range candidate {
			candidate[i] = x[i] + alpha*p[i]
		}
		ax, err := a.Matvec(candidate)
		if err != nil {
			return result("nonfinite")
		}
		res := make([]float64, n)
		for i := range res {
			res[i] = b[i] - ax[i]
		}
		length := norm(res)
		if !isFinite(length) {
			return result("nonfinite")
		}
		x, r = candidate, res
		history = append(history, length)
		if o.Capture {
			frames = append(frames, append([]float64{}, x...))
		}
		if length <= threshold {
			return result("converged")
		}
		for i := range z {
			z[i] = r[i] / diag[i]
		}
		next := dot(r, z)
		if !isFinite(next) {
			return result("nonfinite")
		}
		if next <= 0 {
			return result("breakdown")
		}
		beta := next / rho
		for i := range p {
			p[i] = z[i] + beta*p[i]
		}
		rho = next
	}
	return result("iteration_limit")
}

// ReverseCuthillMcKee returns new-to-old indices for the undirected stored pattern.
func (a *CSRMatrix) ReverseCuthillMcKee() ([]int, error) {
	if a == nil || a.rows != a.cols {
		return nil, errors.New("RCM requires square matrix")
	}
	n := a.rows
	graph := make([]map[int]bool, n)
	for i := range graph {
		graph[i] = map[int]bool{}
	}
	for i := 0; i < n; i++ {
		for k := a.rp[i]; k < a.rp[i+1]; k++ {
			j := a.ci[k]
			if i != j {
				graph[i][j] = true
				graph[j][i] = true
			}
		}
	}
	less := func(i, j int) bool {
		if len(graph[i]) == len(graph[j]) {
			return i < j
		}
		return len(graph[i]) < len(graph[j])
	}
	starts := make([]int, n)
	for i := range starts {
		starts[i] = i
	}
	sort.Slice(starts, func(i, j int) bool { return less(starts[i], starts[j]) })
	seen := make([]bool, n)
	order := []int{}
	for _, start := range starts {
		if seen[start] {
			continue
		}
		q := []int{start}
		seen[start] = true
		for h := 0; h < len(q); h++ {
			next := []int{}
			for j := range graph[q[h]] {
				if !seen[j] {
					next = append(next, j)
				}
			}
			sort.Slice(next, func(i, j int) bool { return less(next[i], next[j]) })
			for _, j := range next {
				seen[j] = true
				q = append(q, j)
			}
		}
		order = append(order, q...)
	}
	for i, j := 0, len(order)-1; i < j; i, j = i+1, j-1 {
		order[i], order[j] = order[j], order[i]
	}
	return order, nil
}
func PermuteVector(p []int, x []float64, inverse bool) ([]float64, error) {
	n := len(x)
	if len(p) != n {
		return nil, errors.New("invalid permutation length")
	}
	seen := make([]bool, n)
	out := make([]float64, n)
	for i, j := range p {
		if j < 0 || j >= n || seen[j] || !isFinite(x[i]) {
			return nil, errors.New("invalid permutation or vector")
		}
		seen[j] = true
	}
	for i, j := range p {
		if inverse {
			out[j] = x[i]
		} else {
			out[i] = x[j]
		}
	}
	return out, nil
}
func (a *CSRMatrix) PermuteSymmetric(p []int) (*CSRMatrix, error) {
	if a == nil || a.rows != a.cols {
		return nil, errors.New("permutation requires square matrix")
	}
	ids := make([]float64, a.rows)
	for i := range ids {
		ids[i] = float64(i)
	}
	inv, e := PermuteVector(p, ids, true)
	if e != nil {
		return nil, e
	}
	rp := []int{0}
	ci := []int{}
	v := []float64{}
	type entry struct {
		j int
		v float64
	}
	for _, i := range p {
		entries := []entry{}
		for k := a.rp[i]; k < a.rp[i+1]; k++ {
			entries = append(entries, entry{int(inv[a.ci[k]]), a.values[k]})
		}
		sort.Slice(entries, func(i, j int) bool { return entries[i].j < entries[j].j })
		for _, x := range entries {
			ci = append(ci, x.j)
			v = append(v, x.v)
		}
		rp = append(rp, len(v))
	}
	return NewCSR(a.rows, a.cols, rp, ci, v)
}
