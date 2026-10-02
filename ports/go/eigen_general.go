package matrix

import (
	"errors"
	"math"
	"sort"
)

// GeneralEigenDecomposition contains lexicographically sorted complex eigenvalues
// and their corresponding unit right eigenvector columns in split real/imaginary form.
// Defective matrices can produce dependent columns; orthogonality is not promised.
type GeneralEigenDecomposition struct {
	EigenvaluesReal, EigenvaluesImag   []float64
	EigenvectorsReal, EigenvectorsImag *Matrix
}

// EigenGeneral solves any finite real square matrix, with 1000 QR iterations per root.
func (m *Matrix) EigenGeneral() (*GeneralEigenDecomposition, error) { return m.EigenGeneralWith(1000) }

// EigenGeneralWith sets the per-root QR iteration limit (1 through 100000).
// Reduction uses scaled Hessenberg double-shift QR and machine-epsilon deflation.
// The input is unchanged. Nonconvergence and nonfinite results return errors.
func (m *Matrix) EigenGeneralWith(maxIterations int) (*GeneralEigenDecomposition, error) {
	if m == nil || m.rows != m.cols {
		return nil, errors.New("eigendecomposition requires a square matrix")
	}
	if maxIterations < 1 || maxIterations > 100000 {
		return nil, errors.New("iteration limit must be between 1 and 100000")
	}
	n := m.rows
	scale := 0.0
	diagonal := true
	upper, lower := true, true
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			value := m.values[i*n+j]
			if !isFinite(value) {
				return nil, errors.New("eigendecomposition requires finite entries")
			}
			scale = math.Max(scale, math.Abs(value))
			if i != j && value != 0 {
				diagonal = false
			}
			if i > j && value != 0 {
				upper = false
			}
			if i < j && value != 0 {
				lower = false
			}
		}
	}
	H, V := make([][]float64, n), make([][]float64, n)
	d, e := make([]float64, n), make([]float64, n)
	for i := 0; i < n; i++ {
		H[i] = make([]float64, n)
		V[i] = make([]float64, n)
		V[i][i] = 1
		copy(H[i], m.values[i*n:(i+1)*n])
		d[i] = H[i][i]
	}
	reverse := lower && !upper
	balance := make([]int, n)
	if diagonal {
		scale = 1
	} else {
		if reverse {
			for i := 0; i < n; i++ {
				for j := 0; j < n; j++ {
					H[i][j] = m.values[(n-1-i)*n+n-1-j]
				}
			}
		}
		if !upper && !lower {
			balance = generalBalance(H)
		}
		scale = 0
		for i := range H {
			for _, v := range H[i] {
				scale = math.Max(scale, math.Abs(v))
			}
		}
		exponent := (math.Float64bits(scale) >> 52) & 0x7ff
		if exponent != 0 {
			scale = math.Float64frombits(exponent << 52)
		}
		for i := range H {
			for j := range H[i] {
				H[i][j] /= scale
			}
		}
		generalHessenberg(H, V, make([]float64, n), n)
		if err := generalQR(H, V, d, e, maxIterations); err != nil {
			return nil, err
		}
	}
	vr, _ := Zeros(n, n)
	vi, _ := Zeros(n, n)
	order := make([]int, n)
	imagSigns := append([]float64(nil), e...)
	for i := 0; i < n; i++ {
		if upper || lower {
			original := i
			if reverse {
				original = n - 1 - i
			}
			d[i] = m.values[original*n+original]
		} else {
			d[i] *= scale
		}
		e[i] *= scale
		if !isFinite(d[i]) || !isFinite(e[i]) {
			return nil, errors.New("eigenvalue is outside the finite float64 range")
		}
		order[i] = i
	}
	sort.SliceStable(order, func(i, j int) bool { a, b := order[i], order[j]; return d[a] < d[b] || d[a] == d[b] && e[a] < e[b] })
	dr, di := make([]float64, n), make([]float64, n)
	for col, original := range order {
		realColumn, imagColumn, imagSign := original, -1, 0.0
		if imagSigns[original] > 0 {
			imagColumn = original + 1
			imagSign = 1
		} else if imagSigns[original] < 0 {
			realColumn = original - 1
			imagColumn = original
			imagSign = -1
		}
		length := 0.0
		for row := 0; row < n; row++ {
			inputRow := row
			if reverse {
				inputRow = n - 1 - row
			}
			r := V[inputRow][realColumn]
			im := 0.0
			if imagColumn >= 0 {
				im = imagSign * V[inputRow][imagColumn]
			}
			if !isFinite(r) || !isFinite(im) {
				return nil, errors.New("nonfinite eigenvector")
			}
			vr.values[row*n+col] = r
			vi.values[row*n+col] = im

		}
		generalUnbalance(vr, vi, col, balance)
		for row := 0; row < n; row++ {
			length = math.Hypot(length, math.Hypot(vr.values[row*n+col], vi.values[row*n+col]))
		}
		if length == 0 || !isFinite(length) {
			return nil, errors.New("invalid eigenvector norm")
		}
		for row := 0; row < n; row++ {
			vr.values[row*n+col] /= length
			vi.values[row*n+col] /= length
		}
		dr[col], di[col] = d[original], e[original]
	}
	return &GeneralEigenDecomposition{dr, di, vr, vi}, nil
}

// Hessenberg reduction and real double-shift QR follow the public-domain JAMA
// implementation of the EISPACK orthes/hqr2 algorithms.
// Reference: https://math.nist.gov/javanumerics/jama/
func generalHessenberg(H, V [][]float64, ort []float64, n int) {

	low := 0
	high := n - 1

	for m := low + 1; m <= high-1; m++ {

		var scale float64
		scale = 0.0
		for i := m; i <= high; i++ {
			scale = scale + math.Abs(H[i][m-1])
		}
		if scale != 0.0 {

			var h float64
			h = 0.0
			for i := high; i >= m; i-- {
				ort[i] = H[i][m-1] / scale
				h += ort[i] * ort[i]
			}
			var g float64
			g = math.Sqrt(h)
			if ort[m] > 0 {
				g = -g
			}
			h = h - ort[m]*g
			ort[m] = ort[m] - g

			for j := m; j < n; j++ {
				var f float64
				f = 0.0
				for i := high; i >= m; i-- {
					f += ort[i] * H[i][j]
				}
				f = f / h
				for i := m; i <= high; i++ {
					H[i][j] -= f * ort[i]
				}
			}

			for i := 0; i <= high; i++ {
				var f float64
				f = 0.0
				for j := high; j >= m; j-- {
					f += ort[j] * H[i][j]
				}
				f = f / h
				for j := m; j <= high; j++ {
					H[i][j] -= f * ort[j]
				}
			}
			ort[m] = scale * ort[m]
			H[m][m-1] = scale * g
		}
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			V[i][j] = 0.0
			if i == j {
				V[i][j] = 1.0
			}
		}
	}

	for m := high - 1; m >= low+1; m-- {
		if H[m][m-1] != 0.0 {
			for i := m + 1; i <= high; i++ {
				ort[i] = H[i][m-1]
			}
			for j := m; j <= high; j++ {
				var g float64
				g = 0.0
				for i := m; i <= high; i++ {
					g += ort[i] * V[i][j]
				}

				g = (g / ort[m]) / H[m][m-1]
				for i := m; i <= high; i++ {
					V[i][j] += g * ort[i]
				}
			}
		}
	}

}

func generalQR(H, V [][]float64, d, e []float64, maxIterations int) error {

	cdivr, cdivi := 0.0, 0.0
	cdiv := func(xr, xi, yr, yi float64) {
		if math.Abs(yr) > math.Abs(yi) {
			r := yi / yr
			den := yr + r*yi
			cdivr = (xr + r*xi) / den
			cdivi = (xi - r*xr) / den
		} else {
			r := yr / yi
			den := yi + r*yr
			cdivr = (r*xr + xi) / den
			cdivi = (r*xi - xr) / den
		}
	}

	nn := len(H)
	n := nn - 1
	low := 0
	high := nn - 1
	var eps float64
	eps = math.Pow(2.0, -52.0)
	var exshift float64
	exshift = 0.0
	var p, q, r, s, z, t, w, x, y float64
	p = 0
	q = 0
	r = 0
	s = 0
	z = 0

	var norm float64
	norm = 0.0
	for i := 0; i < nn; i++ {
		if i < low || i > high {
			d[i] = H[i][i]
			e[i] = 0.0
		}
		for j := max(i-1, 0); j < nn; j++ {
			norm = norm + math.Abs(H[i][j])
		}
	}

	iter := 0
	for n >= low {

		l := n
		for l > low {
			s = math.Abs(H[l-1][l-1]) + math.Abs(H[l][l])
			if s == 0.0 {
				s = norm
			}
			if math.Abs(H[l][l-1]) <= eps*s {
				break
			}
			l--
		}

		if l == n {
			H[n][n] = H[n][n] + exshift
			d[n] = H[n][n]
			e[n] = 0.0
			n--
			iter = 0

		} else if l == n-1 {
			w = H[n][n-1] * H[n-1][n]
			p = (H[n-1][n-1] - H[n][n]) / 2.0
			q = p*p + w
			z = math.Sqrt(math.Abs(q))
			H[n][n] = H[n][n] + exshift
			H[n-1][n-1] = H[n-1][n-1] + exshift
			x = H[n][n]

			if q >= 0 {
				if p >= 0 {
					z = p + z
				} else {
					z = p - z
				}
				d[n-1] = x + z
				d[n] = d[n-1]
				if z != 0.0 {
					d[n] = x - w/z
				}
				e[n-1] = 0.0
				e[n] = 0.0
				x = H[n][n-1]
				s = math.Abs(x) + math.Abs(z)
				p = x / s
				q = z / s
				r = math.Sqrt(p*p + q*q)
				p = p / r
				q = q / r

				for j := n - 1; j < nn; j++ {
					z = H[n-1][j]
					H[n-1][j] = q*z + p*H[n][j]
					H[n][j] = q*H[n][j] - p*z
				}

				for i := 0; i <= n; i++ {
					z = H[i][n-1]
					H[i][n-1] = q*z + p*H[i][n]
					H[i][n] = q*H[i][n] - p*z
				}

				for i := low; i <= high; i++ {
					z = V[i][n-1]
					V[i][n-1] = q*z + p*V[i][n]
					V[i][n] = q*V[i][n] - p*z
				}

			} else {
				d[n-1] = x + p
				d[n] = x + p
				e[n-1] = z
				e[n] = -z
			}
			n = n - 2
			iter = 0

		} else {

			x = H[n][n]
			y = 0.0
			w = 0.0
			if l < n {
				y = H[n-1][n-1]
				w = H[n][n-1] * H[n-1][n]
			}

			if iter == 10 {
				exshift += x
				for i := low; i <= n; i++ {
					H[i][i] -= x
				}
				s = math.Abs(H[n][n-1]) + math.Abs(H[n-1][n-2])
				x = 0.75 * s
				y = x
				w = -0.4375 * s * s
			}

			if iter == 30 {
				s = (y - x) / 2.0
				s = s*s + w
				if s > 0 {
					s = math.Sqrt(s)
					if y < x {
						s = -s
					}
					s = x - w/((y-x)/2.0+s)
					for i := low; i <= n; i++ {
						H[i][i] -= s
					}
					exshift += s
					x = 0.964
					y = x
					w = x
				}
			}

			iter = iter + 1
			if iter > maxIterations {
				return errors.New("general eigendecomposition did not converge within the iteration limit")
			}

			m := n - 2
			for m >= l {
				z = H[m][m]
				r = x - z
				s = y - z
				p = (r*s-w)/H[m+1][m] + H[m][m+1]
				q = H[m+1][m+1] - z - r - s
				r = H[m+2][m+1]
				s = math.Abs(p) + math.Abs(q) + math.Abs(r)
				p = p / s
				q = q / s
				r = r / s
				if m == l {
					break
				}
				if math.Abs(H[m][m-1])*(math.Abs(q)+math.Abs(r)) <
					eps*(math.Abs(p)*(math.Abs(H[m-1][m-1])+math.Abs(z)+
						math.Abs(H[m+1][m+1]))) {
					break
				}
				m--
			}

			for i := m + 2; i <= n; i++ {
				H[i][i-2] = 0.0
				if i > m+2 {
					H[i][i-3] = 0.0
				}
			}

			for k := m; k <= n-1; k++ {
				notlast := (k != n-1)
				if k != m {
					p = H[k][k-1]
					q = H[k+1][k-1]
					r = 0.0
					if notlast {
						r = H[k+2][k-1]
					}
					x = math.Abs(p) + math.Abs(q) + math.Abs(r)
					if x == 0.0 {
						continue
					}
					p = p / x
					q = q / x
					r = r / x
				}

				s = math.Sqrt(p*p + q*q + r*r)
				if p < 0 {
					s = -s
				}
				if s != 0 {
					if k != m {
						H[k][k-1] = -s * x
					} else if l != m {
						H[k][k-1] = -H[k][k-1]
					}
					p = p + s
					x = p / s
					y = q / s
					z = r / s
					q = q / p
					r = r / p

					for j := k; j < nn; j++ {
						p = H[k][j] + q*H[k+1][j]
						if notlast {
							p = p + r*H[k+2][j]
							H[k+2][j] = H[k+2][j] - p*z
						}
						H[k][j] = H[k][j] - p*x
						H[k+1][j] = H[k+1][j] - p*y
					}

					for i := 0; i <= min(n, k+3); i++ {
						p = x*H[i][k] + y*H[i][k+1]
						if notlast {
							p = p + z*H[i][k+2]
							H[i][k+2] = H[i][k+2] - p*r
						}
						H[i][k] = H[i][k] - p
						H[i][k+1] = H[i][k+1] - p*q
					}

					for i := low; i <= high; i++ {
						p = x*V[i][k] + y*V[i][k+1]
						if notlast {
							p = p + z*V[i][k+2]
							V[i][k+2] = V[i][k+2] - p*r
						}
						V[i][k] = V[i][k] - p
						V[i][k+1] = V[i][k+1] - p*q
					}
				}
			}
		}
	}

	if norm == 0.0 {
		return nil
	}

	for n = nn - 1; n >= 0; n-- {
		p = d[n]
		q = e[n]

		if q == 0 {
			l := n
			H[n][n] = 1.0
			for i := n - 1; i >= 0; i-- {
				w = H[i][i] - p
				r = 0.0
				for j := l; j <= n; j++ {
					r = r + H[i][j]*H[j][n]
				}
				if e[i] < 0.0 {
					z = w
					s = r
				} else {
					l = i
					if e[i] == 0.0 {
						if w != 0.0 {
							H[i][n] = -r / w
						} else {
							H[i][n] = -r / (eps * norm)
						}

					} else {
						x = H[i][i+1]
						y = H[i+1][i]
						q = (d[i]-p)*(d[i]-p) + e[i]*e[i]
						t = (x*s - z*r) / q
						H[i][n] = t
						if math.Abs(x) > math.Abs(z) {
							H[i+1][n] = (-r - w*t) / x
						} else {
							H[i+1][n] = (-s - y*t) / z
						}
					}

					t = math.Abs(H[i][n])
					if (eps*t)*t > 1 {
						for j := i; j <= n; j++ {
							H[j][n] = H[j][n] / t
						}
					}
				}
			}

		} else if q < 0 {
			l := n - 1

			if math.Abs(H[n][n-1]) > math.Abs(H[n-1][n]) {
				H[n-1][n-1] = q / H[n][n-1]
				H[n-1][n] = -(H[n][n] - p) / H[n][n-1]
			} else {
				cdiv(0.0, -H[n-1][n], H[n-1][n-1]-p, q)
				H[n-1][n-1] = cdivr
				H[n-1][n] = cdivi
			}
			H[n][n-1] = 0.0
			H[n][n] = 1.0
			for i := n - 2; i >= 0; i-- {
				var ra, sa, vr, vi float64
				ra = 0.0
				sa = 0.0
				for j := l; j <= n; j++ {
					ra = ra + H[i][j]*H[j][n-1]
					sa = sa + H[i][j]*H[j][n]
				}
				w = H[i][i] - p

				if e[i] < 0.0 {
					z = w
					r = ra
					s = sa
				} else {
					l = i
					if e[i] == 0 {
						cdiv(-ra, -sa, w, q)
						H[i][n-1] = cdivr
						H[i][n] = cdivi
					} else {

						x = H[i][i+1]
						y = H[i+1][i]
						vr = (d[i]-p)*(d[i]-p) + e[i]*e[i] - q*q
						vi = (d[i] - p) * 2.0 * q
						if vr == 0.0 && vi == 0.0 {
							vr = eps * norm * (math.Abs(w) + math.Abs(q) +
								math.Abs(x) + math.Abs(y) + math.Abs(z))
						}
						cdiv(x*r-z*ra+q*sa, x*s-z*sa-q*ra, vr, vi)
						H[i][n-1] = cdivr
						H[i][n] = cdivi
						if math.Abs(x) > (math.Abs(z) + math.Abs(q)) {
							H[i+1][n-1] = (-ra - w*H[i][n-1] + q*H[i][n]) / x
							H[i+1][n] = (-sa - w*H[i][n] - q*H[i][n-1]) / x
						} else {
							cdiv(-r-y*H[i][n-1], -s-y*H[i][n], z, q)
							H[i+1][n-1] = cdivr
							H[i+1][n] = cdivi
						}
					}

					t = max(math.Abs(H[i][n-1]), math.Abs(H[i][n]))
					if (eps*t)*t > 1 {
						for j := i; j <= n; j++ {
							H[j][n-1] = H[j][n-1] / t
							H[j][n] = H[j][n] / t
						}
					}
				}
			}
		}
	}

	for i := 0; i < nn; i++ {
		if i < low || i > high {
			for j := i; j < nn; j++ {
				V[i][j] = H[i][j]
			}
		}
	}

	for j := nn - 1; j >= low; j-- {
		for i := low; i <= high; i++ {
			z = 0.0
			for k := low; k <= min(j, high); k++ {
				z = z + V[i][k]*H[k][j]
			}
			V[i][j] = z
		}
	}

	return nil
}

// Balance B=D^-1*A*D without materializing potentially overflowing diagonal D.
func generalBalance(h [][]float64) []int {
	n := len(h)
	exponents := make([]int, n)
	logNorm := func(index int, row bool) float64 {
		largest := 0.0
		for j := 0; j < n; j++ {
			if j != index {
				v := h[index][j]
				if !row {
					v = h[j][index]
				}
				largest = math.Max(largest, math.Abs(v))
			}
		}
		if largest == 0 {
			return math.Inf(-1)
		}
		sum := 0.0
		for j := 0; j < n; j++ {
			if j != index {
				v := h[index][j]
				if !row {
					v = h[j][index]
				}
				sum += math.Abs(v) / largest
			}
		}
		return math.Log2(largest) + math.Log2(sum)
	}
	logAdd := func(a, b float64) float64 {
		top := math.Max(a, b)
		return top + math.Log2(math.Exp2(a-top)+math.Exp2(b-top))
	}
	for sweep := 0; sweep < 64; sweep++ {
		changed := false
		for i := 0; i < n; i++ {
			r, c := logNorm(i, true), logNorm(i, false)
			if math.IsInf(r, -1) || math.IsInf(c, -1) {
				continue
			}
			shift := int(math.Floor((r-c)/2 + .5))
			shift = max(-512, min(512, shift))
			if shift == 0 || logAdd(r-float64(shift), c+float64(shift)) >= logAdd(r, c)+math.Log2(.95) {
				continue
			}
			safe := true
			for j := 0; j < n && safe; j++ {
				if i != j {
					for _, pair := range []struct {
						value float64
						power int
					}{{h[i][j], -shift}, {h[j][i], shift}} {
						scaled := math.Ldexp(pair.value, pair.power)
						if !isFinite(scaled) || (pair.value != 0 && (scaled == 0 || math.Ldexp(scaled, -pair.power) != pair.value)) {
							safe = false
						}
					}
				}
			}
			if !safe {
				continue
			}
			for j := 0; j < n; j++ {
				if i != j {
					h[i][j] = math.Ldexp(h[i][j], -shift)
					h[j][i] = math.Ldexp(h[j][i], shift)
				}
			}
			exponents[i] += shift
			changed = true
		}
		if !changed {
			break
		}
	}
	return exponents
}

func generalUnbalance(vr, vi *Matrix, col int, exponents []int) {
	n := vr.rows
	common := math.MinInt
	for row := 0; row < n; row++ {
		v := math.Max(math.Abs(vr.values[row*n+col]), math.Abs(vi.values[row*n+col]))
		if v != 0 {
			_, e := math.Frexp(v)
			common = max(common, e-1+exponents[row])
		}
	}
	if common == math.MinInt {
		return
	}
	for row := 0; row < n; row++ {
		vr.values[row*n+col] = math.Ldexp(vr.values[row*n+col], exponents[row]-common)
		vi.values[row*n+col] = math.Ldexp(vi.values[row*n+col], exponents[row]-common)
	}
}
