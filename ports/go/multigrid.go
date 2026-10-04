package matrix

import "errors"

// SymmetricPreconditioner must be fixed, SPD, and preserve the input vector.
type SymmetricPreconditioner interface {
	Size() int
	Apply([]float64) ([]float64, error)
}

// GeometricMultigrid is a symmetric V-cycle for the unit Dirichlet Laplacian.
type GeometricMultigrid struct{ width int }

func NewGeometricMultigrid(width int) (*GeometricMultigrid, error) {
	if width < 1 || width > 255 || width&(width+1) != 0 {
		return nil, errors.New("grid width must be 2^k-1 in 1..255")
	}
	return &GeometricMultigrid{width}, nil
}
func (m *GeometricMultigrid) Width() int { return m.width }
func (m *GeometricMultigrid) Size() int  { return m.width * m.width }
func (m *GeometricMultigrid) Levels() int {
	k := 0
	for w := m.width; w > 0; w /= 2 {
		k++
	}
	return k
}
func (m *GeometricMultigrid) Matrix() (*CSRMatrix, error) {
	if m == nil || m.width == 0 {
		return nil, errors.New("uninitialized multigrid")
	}
	w, n := m.width, m.Size()
	rp, ci, v := []int{0}, []int{}, []float64{}
	for i := 0; i < n; i++ {
		add := func(j int) {
			ci = append(ci, j)
			if j == i {
				v = append(v, 4)
			} else {
				v = append(v, -1)
			}
		}
		if i >= w {
			add(i - w)
		}
		if i%w > 0 {
			add(i - 1)
		}
		add(i)
		if i%w+1 < w {
			add(i + 1)
		}
		if i+w < n {
			add(i + w)
		}
		rp = append(rp, len(v))
	}
	return NewCSR(n, n, rp, ci, v)
}
func mgMultiply(w int, x []float64) []float64 {
	r := make([]float64, len(x))
	for i, v := range x {
		r[i] = 4 * v
		if i%w > 0 {
			r[i] -= x[i-1]
		}
		if i%w+1 < w {
			r[i] -= x[i+1]
		}
		if i >= w {
			r[i] -= x[i-w]
		}
		if i+w < len(x) {
			r[i] -= x[i+w]
		}
	}
	return r
}
func mgCycle(w int, b []float64) ([]float64, error) {
	if w == 1 {
		return []float64{b[0] / 4}, nil
	}
	x := make([]float64, len(b))
	smooth := func() {
		for k := 0; k < 2; k++ {
			ax := mgMultiply(w, x)
			for i := range x {
				x[i] += (b[i] - ax[i]) / 6
			}
		}
	}
	smooth()
	r := mgMultiply(w, x)
	for i := range r {
		r[i] = b[i] - r[i]
	}
	c := w / 2
	bc := make([]float64, c*c)
	weight := func(d int) float64 {
		if d == 0 {
			return 1
		}
		return .5
	}
	for y := 0; y < c; y++ {
		for j := 0; j < c; j++ {
			for dy := -1; dy <= 1; dy++ {
				for dx := -1; dx <= 1; dx++ {
					bc[y*c+j] += weight(dy) * weight(dx) * r[(2*y+1+dy)*w+2*j+1+dx]
				}
			}
		}
	}
	ec, e := mgCycle(c, bc)
	if e != nil {
		return nil, e
	}
	for y := 0; y < c; y++ {
		for j := 0; j < c; j++ {
			for dy := -1; dy <= 1; dy++ {
				for dx := -1; dx <= 1; dx++ {
					x[(2*y+1+dy)*w+2*j+1+dx] += weight(dy) * weight(dx) * ec[y*c+j]
				}
			}
		}
	}
	smooth()
	for _, v := range x {
		if !isFinite(v) {
			return nil, errors.New("nonfinite multigrid cycle")
		}
	}
	return x, nil
}
func (m *GeometricMultigrid) Apply(b []float64) ([]float64, error) {
	if m == nil || m.width == 0 || len(b) != m.Size() {
		return nil, errors.New("invalid multigrid RHS")
	}
	for _, v := range b {
		if !isFinite(v) {
			return nil, errors.New("nonfinite multigrid RHS")
		}
	}
	return mgCycle(m.width, b)
}
