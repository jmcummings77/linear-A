package matrix

import "errors"

// ILU0 owns a zero-fill factorization with no pivoting or diagonal shifts.
type ILU0 struct {
	a        *CSRMatrix
	diagonal []int
}

func NewILU0(a *CSRMatrix) (*ILU0, error) {
	if a == nil || a.rows != a.cols {
		return nil, errors.New("ILU0 requires square matrix")
	}
	copyA, err := NewCSR(a.rows, a.cols, a.rp, a.ci, a.values)
	if err != nil {
		return nil, err
	}
	a = copyA
	d := make([]int, a.rows)
	for i := range d {
		d[i] = a.find(i, i)
		if d[i] == a.rp[i+1] || a.ci[d[i]] != i {
			return nil, errors.New("ILU0 requires stored diagonal")
		}
	}
	for i := 0; i < a.rows; i++ {
		for p := a.rp[i]; p < d[i]; p++ {
			j := a.ci[p]
			a.values[p] /= a.values[d[j]]
			if !isFinite(a.values[p]) {
				return nil, errors.New("nonfinite ILU0 factor")
			}
			for q := d[j] + 1; q < a.rp[j+1]; q++ {
				k := a.find(i, a.ci[q])
				if k < a.rp[i+1] && a.ci[k] == a.ci[q] {
					a.values[k] -= a.values[p] * a.values[q]
					if !isFinite(a.values[k]) {
						return nil, errors.New("nonfinite ILU0 factor")
					}
				}
			}
		}
		if a.values[d[i]] == 0 {
			return nil, errors.New("zero ILU0 pivot")
		}
	}
	return &ILU0{a, d}, nil
}
func (f *ILU0) Size() int { return f.a.rows }
func (f *ILU0) NNZ() int  { return len(f.a.values) }
func (f *ILU0) Apply(b []float64) ([]float64, error) {
	if f == nil || f.a == nil || len(b) != f.a.rows {
		return nil, errors.New("invalid ILU0 vector")
	}
	x := append([]float64{}, b...)
	for _, v := range x {
		if !isFinite(v) {
			return nil, errors.New("nonfinite ILU0 vector")
		}
	}
	a, d := f.a, f.diagonal
	for i := 0; i < a.rows; i++ {
		for p := a.rp[i]; p < d[i]; p++ {
			x[i] -= a.values[p] * x[a.ci[p]]
		}
		if !isFinite(x[i]) {
			return nil, errors.New("nonfinite ILU0 solve")
		}
	}
	for i := a.rows - 1; i >= 0; i-- {
		for p := d[i] + 1; p < a.rp[i+1]; p++ {
			x[i] -= a.values[p] * x[a.ci[p]]
		}
		x[i] /= a.values[d[i]]
		if !isFinite(x[i]) {
			return nil, errors.New("nonfinite ILU0 solve")
		}
	}
	return x, nil
}
