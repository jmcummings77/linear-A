package main

import (
	"encoding/json"
	"errors"
	matrix "github.com/jmcummings77/linear-A/ports/go"
	"io"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

func sparseRun(args []string) error {
	expected := 11
	if len(args) > 1 && args[1] == "gmres" {
		expected = 12
	}
	if len(args) != expected {
		return errors.New("invalid sparse protocol")
	}
	op := args[1]
	ints := []int{}
	for _, index := range []int{2, 3, 4, 5, 8, 9, 10} {
		v, e := strconv.Atoi(args[index])
		if e != nil || v < 0 {
			return errors.New("invalid sparse integer")
		}
		ints = append(ints, v)
	}
	rows, cols, nnz, iterations, limit, jacobi, capture := ints[0], ints[1], ints[2], ints[3], ints[4], ints[5], ints[6]
	rtol, e := strconv.ParseFloat(args[6], 64)
	if e != nil {
		return e
	}
	atol, e := strconv.ParseFloat(args[7], 64)
	if e != nil {
		return e
	}
	if (op != "mg_setup" && op != "mg_matrix" && op != "mg_apply" && op != "spmv" && op != "dense" && op != "cg" && op != "gmres" && op != "ic0_factor" && op != "ic0_apply" && op != "ilu_setup" && op != "ilu_apply" && op != "rcm" && op != "amd" && op != "permute" && op != "permutation_check" && op != "rcm_solve" && op != "ilu_solve" && op != "chol_symbolic" && op != "chol_factor" && op != "chol_solve" && op != "chol_total" && op != "chol_rcm_total" && op != "chol_amd_total") || jacobi > 5 || (op != "cg" && jacobi > 3) || (op != "gmres" && op != "cg" && jacobi > 1) || capture > 1 {
		return errors.New("invalid sparse operation/options")
	}
	raw, e := io.ReadAll(os.Stdin)
	if e != nil {
		return e
	}
	tokens := strings.Fields(string(raw))
	count := cols
	if op == "cg" || op == "gmres" {
		count = rows
	}
	if len(tokens) != rows+1+2*nnz+count {
		return errors.New("incorrect sparse input count")
	}
	offset := 0
	takeInt := func(n int) ([]int, error) {
		a := make([]int, n)
		for i := range a {
			v, e := strconv.Atoi(tokens[offset])
			if e != nil {
				return nil, e
			}
			offset++
			a[i] = v
		}
		return a, nil
	}
	take := func(n int) ([]float64, error) {
		a := make([]float64, n)
		for i := range a {
			v, e := strconv.ParseFloat(tokens[offset], 64)
			if e != nil {
				return nil, e
			}
			offset++
			a[i] = v
		}
		return a, nil
	}
	rp, e := takeInt(rows + 1)
	if e != nil {
		return e
	}
	ci, e := takeInt(nnz)
	if e != nil {
		return e
	}
	v, e := take(nnz)
	if e != nil {
		return e
	}
	b, e := take(count)
	if e != nil {
		return e
	}
	a, e := matrix.NewCSR(rows, cols, rp, ci, v)
	if e != nil {
		return e
	}
	newMG := func() (*matrix.GeometricMultigrid, error) {
		w := int(math.Sqrt(float64(rows)))
		if rows != cols || w*w != rows {
			return nil, errors.New("multigrid requires a square grid")
		}
		return matrix.NewGeometricMultigrid(w)
	}
	var mg *matrix.GeometricMultigrid
	if op == "mg_apply" || op == "mg_matrix" || (op == "cg" && jacobi == 4) {
		mg, e = newMG()
		if e != nil {
			return e
		}
	}
	var ic *matrix.IC0
	if op == "ic0_apply" || (op == "cg" && jacobi == 2) {
		ic, e = matrix.NewIC0(a)
		if e != nil {
			return e
		}
	}
	var factor *matrix.ILU0
	if op == "ilu_apply" || (op == "gmres" && jacobi == 2) {
		factor, e = matrix.NewILU0(a)
		if e != nil {
			return e
		}
	}
	var dense, right *matrix.Matrix
	if op == "dense" {
		dv := make([]float64, rows*cols)
		for i := 0; i < rows; i++ {
			for p := rp[i]; p < rp[i+1]; p++ {
				dv[i*cols+ci[p]] = v[p]
			}
		}
		dense, e = matrix.New(rows, cols, dv)
		if e != nil {
			return e
		}
		right, e = matrix.New(cols, 1, b)
		if e != nil {
			return e
		}
	}
	var plan *matrix.SparseCholeskySymbolic
	var chol *matrix.SparseCholesky
	if op == "chol_factor" || op == "chol_solve" {
		plan, e = matrix.NewSparseCholeskySymbolic(a)
		if e != nil {
			return e
		}
	}
	if op == "chol_solve" {
		chol, e = plan.Factorize(a)
		if e != nil {
			return e
		}
	}
	compute := func() ([]float64, error) {
		if op == "mg_setup" {
			m, e := newMG()
			if e != nil {
				return nil, e
			}
			return []float64{float64(m.Size()), float64(m.Levels())}, nil
		}
		if op == "mg_apply" {
			return mg.Apply(b)
		}
		if op == "mg_matrix" {
			m, e := mg.Matrix()
			if e != nil {
				return nil, e
			}
			out := []float64{}
			for _, v := range m.RowOffsets() {
				out = append(out, float64(v))
			}
			for _, v := range m.ColumnIndices() {
				out = append(out, float64(v))
			}
			return append(out, m.Values()...), nil
		}
		if op == "ic0_factor" {
			f, e := matrix.NewIC0(a)
			if e != nil {
				return nil, e
			}
			l := f.Lower()
			out := []float64{}
			for _, v := range l.RowOffsets() {
				out = append(out, float64(v))
			}
			for _, v := range l.ColumnIndices() {
				out = append(out, float64(v))
			}
			return append(out, l.Values()...), nil
		}
		if op == "ic0_apply" {
			return ic.Apply(b)
		}
		if op == "chol_symbolic" {
			s, e := matrix.NewSparseCholeskySymbolic(a)
			if e != nil {
				return nil, e
			}
			v := []float64{}
			for _, xs := range [][]int{s.RowOffsets(), s.ColumnIndices(), s.FillSteps()} {
				for _, x := range xs {
					v = append(v, float64(x))
				}
			}
			return v, nil
		}
		if op == "chol_factor" {
			f, e := plan.Factorize(a)
			if e != nil {
				return nil, e
			}
			l := f.Lower()
			v := []float64{}
			for _, xs := range [][]int{l.RowOffsets(), l.ColumnIndices()} {
				for _, x := range xs {
					v = append(v, float64(x))
				}
			}
			return append(v, l.Values()...), nil
		}
		if op == "chol_solve" {
			return chol.Solve(b)
		}
		if op == "chol_total" || op == "chol_rcm_total" || op == "chol_amd_total" {
			q, rhs := a, b
			var p []int
			var e error
			if op == "chol_rcm_total" || op == "chol_amd_total" {
				if op == "chol_amd_total" {
					p, e = a.ApproximateMinimumDegree()
				} else {
					p, e = a.ReverseCuthillMcKee()
				}
				if e != nil {
					return nil, e
				}
				q, e = a.PermuteSymmetric(p)
				if e != nil {
					return nil, e
				}
				rhs, e = matrix.PermuteVector(p, b, false)
				if e != nil {
					return nil, e
				}
			}
			s, e := matrix.NewSparseCholeskySymbolic(q)
			if e != nil {
				return nil, e
			}
			f, e := s.Factorize(q)
			if e != nil {
				return nil, e
			}
			x, e := f.Solve(rhs)
			if e != nil {
				return nil, e
			}
			if op == "chol_rcm_total" || op == "chol_amd_total" {
				return matrix.PermuteVector(p, x, true)
			}
			return x, nil
		}

		if op == "rcm_solve" || op == "ilu_solve" {
			q := a
			rhs := b
			var p []int
			var err error
			if op == "rcm_solve" {
				p, err = a.ReverseCuthillMcKee()
				if err != nil {
					return nil, err
				}
				q, err = a.PermuteSymmetric(p)
				if err != nil {
					return nil, err
				}
				rhs, err = matrix.PermuteVector(p, b, false)
				if err != nil {
					return nil, err
				}
			}
			f, err := matrix.NewILU0(q)
			if err != nil {
				return nil, err
			}
			o := matrix.DefaultGMRESOptions()
			o.Restart = 20
			o.RelativeTolerance = rtol
			o.AbsoluteTolerance = atol
			o.MaxIterations = limit
			o.Preconditioner = f
			r, err := q.GMRES(rhs, o)
			if err != nil {
				return nil, err
			}
			if !r.Converged {
				return nil, errors.New("ordering solve failed")
			}
			x := r.X
			if op == "rcm_solve" {
				x, err = matrix.PermuteVector(p, x, true)
				if err != nil {
					return nil, err
				}
			}
			return append([]float64{float64(r.Iterations)}, x...), nil
		}
		if op == "amd" {
			p, e := a.ApproximateMinimumDegree()
			if e != nil {
				return nil, e
			}
			out := make([]float64, len(p))
			for i, j := range p {
				out[i] = float64(j)
			}
			return out, nil
		}
		if op == "rcm" {
			p, e := a.ReverseCuthillMcKee()
			if e != nil {
				return nil, e
			}
			out := make([]float64, len(p))
			for i, j := range p {
				out[i] = float64(j)
			}
			return out, nil
		}
		if op == "permute" || op == "permutation_check" {
			p := make([]int, len(b))
			for i, x := range b {
				if math.IsNaN(x) || math.IsInf(x, 0) || x < 0 || x >= float64(rows) || x != math.Floor(x) {
					return nil, errors.New("invalid permutation")
				}
				p[i] = int(x)
			}
			q, e := a.PermuteSymmetric(p)
			if e != nil {
				return nil, e
			}
			if op == "permute" {
				out := []float64{}
				for _, x := range q.RowOffsets() {
					out = append(out, float64(x))
				}
				for _, x := range q.ColumnIndices() {
					out = append(out, float64(x))
				}
				return append(out, q.Values()...), nil
			}
			x := make([]float64, rows)
			for i := range x {
				x[i] = float64(i + 1)
			}
			y, e := matrix.PermuteVector(p, x, false)
			if e != nil {
				return nil, e
			}
			back, e := matrix.PermuteVector(p, y, true)
			if e != nil {
				return nil, e
			}
			z, e := q.Matvec(y)
			if e != nil {
				return nil, e
			}
			z, e = matrix.PermuteVector(p, z, true)
			if e != nil {
				return nil, e
			}
			return append(append(y, back...), z...), nil
		}
		if op == "ilu_setup" {
			f, e := matrix.NewILU0(a)
			if e != nil {
				return nil, e
			}
			return []float64{float64(f.Size()), float64(f.NNZ())}, nil
		}
		if op == "ilu_apply" {
			return factor.Apply(b)
		}
		if op == "spmv" {
			return a.Matvec(b)
		}
		if op == "dense" {
			r, e := dense.Multiply(right)
			if e != nil {
				return nil, e
			}
			return r.Values(), nil
		}

		if op == "gmres" {
			restart, e := strconv.Atoi(args[11])
			if e != nil {
				return nil, e
			}
			current := factor
			if jacobi == 3 {
				current, e = matrix.NewILU0(a)
				if e != nil {
					return nil, e
				}
			}
			r, e := a.GMRES(b, matrix.GMRESOptions{Restart: restart, RelativeTolerance: rtol, AbsoluteTolerance: atol, MaxIterations: limit, Jacobi: jacobi == 1, Capture: capture != 0, Preconditioner: current})
			if e != nil {
				return nil, e
			}
			if iterations > 0 {
				if !r.Converged {
					return nil, errors.New("benchmark GMRES did not converge: " + r.Reason)
				}
				return r.X, nil
			}
			reason := map[string]int{"converged": 0, "iteration_limit": 1, "breakdown": 2, "nonfinite": 3, "stagnation": 4}[r.Reason]
			out := []float64{float64(reason), float64(r.Iterations), float64(len(r.Residuals))}
			out = append(out, r.X...)
			out = append(out, r.Residuals...)
			out = append(out, r.EstimatedResiduals...)
			out = append(out, float64(len(r.Restarts)))
			for _, v := range r.Restarts {
				out = append(out, float64(v))
			}
			for _, f := range r.Iterates {
				out = append(out, f...)
			}
			return out, nil
		}
		var currentIC matrix.SymmetricPreconditioner
		if ic != nil {
			currentIC = ic
		}
		if jacobi == 4 {
			currentIC = mg
		}
		if jacobi == 5 {
			currentIC, e = newMG()
			if e != nil {
				return nil, e
			}
		}
		if jacobi == 3 {
			currentIC, e = matrix.NewIC0(a)
			if e != nil {
				return nil, e
			}
		}
		r, e := a.ConjugateGradient(b, matrix.CGOptions{RelativeTolerance: rtol, AbsoluteTolerance: atol, MaxIterations: limit, Jacobi: jacobi == 1, Capture: capture != 0, Preconditioner: currentIC})
		if e != nil {
			return nil, e
		}
		if iterations > 0 {
			if !r.Converged {
				return nil, errors.New("benchmark CG did not converge: " + r.Reason)
			}
			return r.X, nil
		}
		reason := map[string]int{"converged": 0, "iteration_limit": 1, "breakdown": 2, "nonfinite": 3}[r.Reason]
		out := []float64{float64(reason), float64(r.Iterations), float64(len(r.Residuals))}
		out = append(out, r.X...)
		out = append(out, r.Residuals...)
		for _, f := range r.Iterates {
			out = append(out, f...)
		}
		return out, nil
	}
	if iterations == 0 {
		out, e := compute()
		if e != nil {
			return e
		}
		return json.NewEncoder(os.Stdout).Encode(map[string]interface{}{"rows": 1, "cols": len(out), "values": out})
	}
	for i := 0; i < 3; i++ {
		if _, e = compute(); e != nil {
			return e
		}
	}
	checksum := 0.
	start := time.Now()
	for i := 0; i < iterations; i++ {
		out, e := compute()
		if e != nil {
			return e
		}
		for _, x := range out {
			checksum += x
		}
	}
	elapsed := time.Since(start).Nanoseconds()
	return json.NewEncoder(os.Stdout).Encode(map[string]interface{}{"elapsed_ns": elapsed, "iterations": iterations, "checksum": checksum})
}
