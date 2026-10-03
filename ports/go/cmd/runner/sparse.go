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
	if (op != "spmv" && op != "dense" && op != "cg" && op != "gmres" && op != "ilu_setup" && op != "ilu_apply" && op != "rcm" && op != "permute" && op != "permutation_check" && op != "rcm_solve" && op != "ilu_solve") || jacobi > 3 || (op != "gmres" && jacobi > 1) || capture > 1 {
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
	compute := func() ([]float64, error) {
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
		r, e := a.ConjugateGradient(b, matrix.CGOptions{RelativeTolerance: rtol, AbsoluteTolerance: atol, MaxIterations: limit, Jacobi: jacobi != 0, Capture: capture != 0})
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
