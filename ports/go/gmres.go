package matrix

import (
	"errors"
	"math"
)

type GMRESOptions struct {
	Preconditioner                       *ILU0
	Restart                              int
	RelativeTolerance, AbsoluteTolerance float64
	MaxIterations                        int
	Jacobi, Capture                      bool
}

func DefaultGMRESOptions() GMRESOptions {
	return GMRESOptions{Restart: 30, RelativeTolerance: 1e-10, MaxIterations: 1000}
}

type GMRESResult struct {
	X                             []float64
	Converged                     bool
	Iterations                    int
	Reason                        string
	Residuals, EstimatedResiduals []float64
	Iterates                      [][]float64
	Restarts                      []int
}

func (a *CSRMatrix) GMRES(b []float64, o GMRESOptions) (*GMRESResult, error) {
	n := a.rows
	if a.cols != n || len(b) != n {
		return nil, errors.New("GMRES requires square matrix and matching vector")
	}
	for _, v := range b {
		if !isFinite(v) {
			return nil, errors.New("nonfinite RHS")
		}
	}
	if o.Restart < 1 || o.Restart > 1024 || o.MaxIterations < 0 || o.MaxIterations > 100000 || !isFinite(o.RelativeTolerance) || o.RelativeTolerance < 0 || o.RelativeTolerance >= 1 || !isFinite(o.AbsoluteTolerance) || o.AbsoluteTolerance < 0 {
		return nil, errors.New("invalid GMRES options")
	}
	if o.Preconditioner != nil && (o.Jacobi || o.Preconditioner.a == nil || o.Preconditioner.Size() != n) {
		return nil, errors.New("incompatible preconditioner")
	}
	diag := make([]float64, n)
	for i := range diag {
		diag[i] = 1
		if o.Jacobi {
			p := a.find(i, i)
			if p == a.rp[i+1] || a.ci[p] != i || a.values[p] == 0 {
				return nil, errors.New("Jacobi requires nonzero diagonal")
			}
			diag[i] = a.values[p]
		}
	}
	norm := func(v []float64) float64 {
		s := 0.
		for _, x := range v {
			s = math.Hypot(s, x)
		}
		return s
	}
	finite := func(v []float64) bool {
		for _, x := range v {
			if !isFinite(x) {
				return false
			}
		}
		return true
	}
	x := make([]float64, n)
	r := append([]float64{}, b...)
	history := []float64{norm(r)}
	estimates := append([]float64{}, history...)
	frames := [][]float64{}
	restarts := []int{}
	if o.Capture {
		frames = append(frames, append([]float64{}, x...))
	}
	result := func(reason string) (*GMRESResult, error) {
		return &GMRESResult{append([]float64{}, x...), reason == "converged", len(history) - 1, reason, history, estimates, frames, restarts}, nil
	}
	threshold := math.Max(o.AbsoluteTolerance, o.RelativeTolerance*history[0])
	m := o.Restart
	if m > n {
		m = n
	}
	if m > o.MaxIterations {
		m = o.MaxIterations
	}
	if !isFinite(history[0]) {
		return result("nonfinite")
	}
	if history[0] <= threshold {
		return result("converged")
	}
	for len(history)-1 < o.MaxIterations {
		if len(history) > 1 {
			restarts = append(restarts, len(history)-1)
		}
		base := append([]float64{}, x...)
		beta := norm(r)
		basis := make([][]float64, m+1)
		basis[0] = make([]float64, n)
		for i := range r {
			basis[0][i] = r[i] / beta
		}
		h := make([][]float64, m+1)
		for i := range h {
			h[i] = make([]float64, m)
		}
		cs, sn, g := make([]float64, m), make([]float64, m), make([]float64, m+1)
		g[0] = beta
		steps := m
		if steps > o.MaxIterations-(len(history)-1) {
			steps = o.MaxIterations - (len(history) - 1)
		}
		for j := 0; j < steps; j++ {
			z := make([]float64, n)
			for i := range z {
				z[i] = basis[j][i] / diag[i]
			}
			if o.Preconditioner != nil {
				var err error
				z, err = o.Preconditioner.Apply(basis[j])
				if err != nil {
					return result("nonfinite")
				}
			}
			w, err := a.Matvec(z)
			if err != nil {
				return result("nonfinite")
			}
			original := norm(w)
			for pass := 0; pass < 2; pass++ {
				for k := 0; k <= j; k++ {
					dot := 0.
					for i := range w {
						dot += basis[k][i] * w[i]
					}
					h[k][j] += dot
					for i := range w {
						w[i] -= dot * basis[k][i]
					}
				}
			}
			tail := norm(w)
			if !isFinite(original) || !isFinite(tail) {
				return result("nonfinite")
			}
			for k := 0; k <= j; k++ {
				if !isFinite(h[k][j]) {
					return result("nonfinite")
				}
			}
			happy := tail <= 8*2.220446049250313e-16*original
			if !happy {
				h[j+1][j] = tail
				basis[j+1] = make([]float64, n)
				for i := range w {
					basis[j+1][i] = w[i] / tail
				}
			}
			for k := 0; k < j; k++ {
				top := cs[k]*h[k][j] + sn[k]*h[k+1][j]
				h[k+1][j] = -sn[k]*h[k][j] + cs[k]*h[k+1][j]
				h[k][j] = top
			}
			pivot := math.Hypot(h[j][j], h[j+1][j])
			if !isFinite(pivot) {
				return result("nonfinite")
			}
			if pivot == 0 {
				return result("breakdown")
			}
			cs[j] = h[j][j] / pivot
			sn[j] = h[j+1][j] / pivot
			h[j][j] = pivot
			h[j+1][j] = 0
			g[j+1] = -sn[j] * g[j]
			g[j] = cs[j] * g[j]
			y := append([]float64{}, g[:j+1]...)
			for k := j; k >= 0; k-- {
				if h[k][k] == 0 {
					return result("breakdown")
				}
				sum := 0.
				for q := k + 1; q <= j; q++ {
					sum += h[k][q] * y[q]
				}
				y[k] = (y[k] - sum) / h[k][k]
			}
			candidate := make([]float64, n)
			for i := range candidate {
				sum := 0.
				for k := 0; k <= j; k++ {
					sum += basis[k][i] * y[k]
				}
				if o.Preconditioner != nil {
					candidate[i] = sum
				} else {
					candidate[i] = base[i] + sum/diag[i]
				}
			}
			if o.Preconditioner != nil {
				var err error
				candidate, err = o.Preconditioner.Apply(candidate)
				if err != nil {
					return result("nonfinite")
				}
				for i := range candidate {
					candidate[i] += base[i]
				}
			}
			if !finite(y) || !finite(candidate) || !isFinite(g[j+1]) {
				return result("nonfinite")
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
			estimates = append(estimates, math.Abs(g[j+1]))
			if o.Capture {
				frames = append(frames, append([]float64{}, x...))
			}
			if length <= threshold {
				return result("converged")
			}
			if happy {
				return result("breakdown")
			}
		}
		same := true
		for i := range x {
			if x[i] != base[i] {
				same = false
			}
		}
		if same {
			return result("stagnation")
		}
	}
	return result("iteration_limit")
}
