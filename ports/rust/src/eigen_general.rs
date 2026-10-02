use crate::{Matrix, MatrixError};
use std::ops::{Index, IndexMut};

/// Complex eigenpairs, sorted by (real, imaginary), in split real/imaginary form.
/// Columns are unit right eigenvectors; defective matrices may have dependent
/// columns and general eigenvectors are not promised to be orthogonal.
#[derive(Clone, Debug)]
pub struct GeneralEigenDecomposition {
    pub eigenvalues_real: Vec<f64>,
    pub eigenvalues_imag: Vec<f64>,
    pub eigenvectors_real: Matrix,
    pub eigenvectors_imag: Matrix,
}

impl Matrix {
    /// General real-square eigenpairs using at most 1000 QR iterations per root.
    pub fn eigen_general(&self) -> Result<GeneralEigenDecomposition, MatrixError> {
        self.eigen_general_with(1000)
    }
    /// Bounded scaled Hessenberg double-shift QR with machine-epsilon deflation.
    /// The limit is 1..=100000; input is unchanged and numerical failures return errors.
    pub fn eigen_general_with(
        &self,
        max_iterations: usize,
    ) -> Result<GeneralEigenDecomposition, MatrixError> {
        self.require_square()?;
        if !(1..=100000).contains(&max_iterations) {
            return Err(MatrixError::new(
                "iteration limit must be between 1 and 100000",
            ));
        }
        let n = self.rows;
        let mut scale = 0.0_f64;
        let mut diagonal = true;
        let (mut upper, mut lower) = (true, true);
        for i in 0..n {
            for j in 0..n {
                let value = self.values[i * n + j];
                if !value.is_finite() {
                    return Err(MatrixError::new(
                        "eigendecomposition requires finite entries",
                    ));
                }
                scale = scale.max(value.abs());
                if i != j && value != 0.0 {
                    diagonal = false;
                }
                if i > j && value != 0.0 {
                    upper = false;
                }
                if i < j && value != 0.0 {
                    lower = false;
                }
            }
        }
        let mut hmat = WorkMatrix {
            n,
            data: self.values.clone(),
        };
        let mut vmat = WorkMatrix {
            n,
            data: vec![0.0; n * n],
        };
        let mut d = WorkVector(vec![0.0; n]);
        let mut e = WorkVector(vec![0.0; n]);
        for i in 0..n {
            vmat.data[i * n + i] = 1.0;
            d.0[i] = hmat.data[i * n + i];
        }
        let reverse = lower && !upper;
        let mut balance = vec![0; n];
        if diagonal {
            scale = 1.0;
        } else {
            if reverse {
                for i in 0..n {
                    for j in 0..n {
                        hmat.data[i * n + j] = self.values[(n - 1 - i) * n + n - 1 - j];
                    }
                }
            }
            if !upper && !lower {
                balance = general_balance(&mut hmat);
            }
            scale = hmat.data.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
            let exponent = (scale.to_bits() >> 52) & 0x7ff;
            if exponent != 0 {
                scale = f64::from_bits(exponent << 52);
            }
            for value in &mut hmat.data {
                *value /= scale;
            }
            general_hessenberg(
                &mut hmat,
                &mut vmat,
                &mut WorkVector(vec![0.0; n]),
                n as isize,
            );
            general_qr(
                &mut hmat,
                &mut vmat,
                &mut d,
                &mut e,
                max_iterations as isize,
            )?;
        }
        let imag_signs = e.0.clone();
        for i in 0..n {
            if upper || lower {
                let original = if reverse { n - 1 - i } else { i };
                d.0[i] = self.values[original * n + original];
            } else {
                d.0[i] *= scale;
            }
            e.0[i] *= scale;
            if !d.0[i].is_finite() || !e.0[i].is_finite() {
                return Err(MatrixError::new(
                    "eigenvalue is outside the finite float64 range",
                ));
            }
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            d.0[a]
                .partial_cmp(&d.0[b])
                .unwrap()
                .then(e.0[a].partial_cmp(&e.0[b]).unwrap())
        });
        let mut vr = Matrix::zeros(n, n)?;
        let mut vi = Matrix::zeros(n, n)?;
        for (col, &original) in order.iter().enumerate() {
            let (rc, ic, sign) = if imag_signs[original] > 0.0 {
                (original, Some(original + 1), 1.0)
            } else if imag_signs[original] < 0.0 {
                (original - 1, Some(original), -1.0)
            } else {
                (original, None, 0.0)
            };
            let mut length = 0.0_f64;
            for row in 0..n {
                let input_row = if reverse { n - 1 - row } else { row };
                let r = vmat.data[input_row * n + rc];
                let im = ic.map_or(0.0, |i| sign * vmat.data[input_row * n + i]);
                if !r.is_finite() || !im.is_finite() {
                    return Err(MatrixError::new("nonfinite eigenvector"));
                }
                vr.values[row * n + col] = r;
                vi.values[row * n + col] = im;
            }
            general_unbalance(&mut vr, &mut vi, col, &balance);
            for row in 0..n {
                length = length.hypot(vr.values[row * n + col].hypot(vi.values[row * n + col]));
            }
            if length == 0.0 || !length.is_finite() {
                return Err(MatrixError::new("invalid eigenvector norm"));
            }
            for row in 0..n {
                vr.values[row * n + col] /= length;
                vi.values[row * n + col] /= length;
            }
        }
        Ok(GeneralEigenDecomposition {
            eigenvalues_real: order.iter().map(|&i| d.0[i]).collect(),
            eigenvalues_imag: order.iter().map(|&i| e.0[i]).collect(),
            eigenvectors_real: vr,
            eigenvectors_imag: vi,
        })
    }
}

// Signed workspace indices keep the Hessenberg/QR boundary arithmetic explicit.
struct WorkMatrix {
    n: usize,
    data: Vec<f64>,
}
impl Index<(isize, isize)> for WorkMatrix {
    type Output = f64;
    fn index(&self, (i, j): (isize, isize)) -> &f64 {
        &self.data[i as usize * self.n + j as usize]
    }
}
impl IndexMut<(isize, isize)> for WorkMatrix {
    fn index_mut(&mut self, (i, j): (isize, isize)) -> &mut f64 {
        &mut self.data[i as usize * self.n + j as usize]
    }
}
struct WorkVector(Vec<f64>);
impl Index<isize> for WorkVector {
    type Output = f64;
    fn index(&self, i: isize) -> &f64 {
        &self.0[i as usize]
    }
}
impl IndexMut<isize> for WorkVector {
    fn index_mut(&mut self, i: isize) -> &mut f64 {
        &mut self.0[i as usize]
    }
}
fn general_cdiv(xr: f64, xi: f64, yr: f64, yi: f64) -> (f64, f64) {
    if yr.abs() > yi.abs() {
        let r = yi / yr;
        let den = yr + r * yi;
        ((xr + r * xi) / den, (xi - r * xr) / den)
    } else {
        let r = yr / yi;
        let den = yi + r * yr;
        ((r * xr + xi) / den, (r * xi - xr) / den)
    }
}

// Hessenberg reduction and real double-shift QR follow the public-domain JAMA
// implementation of EISPACK orthes/hqr2: https://math.nist.gov/javanumerics/jama/
#[allow(unused_mut)]
fn general_hessenberg(
    hmat: &mut WorkMatrix,
    vmat: &mut WorkMatrix,
    ort: &mut WorkVector,
    n: isize,
) {
    let mut low: isize = 0;
    let mut high: isize = n - 1;
    for m in low + 1..=high - 1 {
        let mut scale: f64 = 0.0;
        for i in m..=high {
            scale = scale + f64::abs(hmat[(i, m - 1)]);
        }
        if scale != 0.0 {
            let mut h: f64 = 0.0;
            for i in (m..=high).rev() {
                ort[i] = hmat[(i, m - 1)] / scale;
                h += ort[i] * ort[i];
            }
            let mut g: f64 = f64::sqrt(h);
            if ort[m] > 0.0 {
                g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;
            for j in m..n {
                let mut f: f64 = 0.0;
                for i in (m..=high).rev() {
                    f += ort[i] * hmat[(i, j)];
                }
                f = f / h;
                for i in m..=high {
                    hmat[(i, j)] = hmat[(i, j)] - (f * ort[i]);
                }
            }
            for i in 0..=high {
                let mut f: f64 = 0.0;
                for j in (m..=high).rev() {
                    f += ort[j] * hmat[(i, j)];
                }
                f = f / h;
                for j in m..=high {
                    hmat[(i, j)] = hmat[(i, j)] - (f * ort[j]);
                }
            }
            ort[m] = scale * ort[m];
            hmat[(m, m - 1)] = scale * g;
        }
    }
    for i in 0..n {
        for j in 0..n {
            vmat[(i, j)] = 0.0;
            if i == j {
                vmat[(i, j)] = 1.0;
            }
        }
    }
    for m in (low + 1..=high - 1).rev() {
        if hmat[(m, m - 1)] != 0.0 {
            for i in m + 1..=high {
                ort[i] = hmat[(i, m - 1)];
            }
            for j in m..=high {
                let mut g: f64 = 0.0;
                for i in m..=high {
                    g += ort[i] * vmat[(i, j)];
                }
                g = (g / ort[m]) / hmat[(m, m - 1)];
                for i in m..=high {
                    vmat[(i, j)] = vmat[(i, j)] + (g * ort[i]);
                }
            }
        }
    }
}

#[allow(unused_mut, unused_assignments)]
fn general_qr(
    hmat: &mut WorkMatrix,
    vmat: &mut WorkMatrix,
    d: &mut WorkVector,
    e: &mut WorkVector,
    max_iterations: isize,
) -> Result<(), MatrixError> {
    let (mut cdivr, mut cdivi) = (0.0, 0.0);
    let mut nn: isize = hmat.n as isize;
    let mut n: isize = nn - 1;
    let mut low: isize = 0;
    let mut high: isize = nn - 1;
    let mut eps: f64 = f64::EPSILON;
    let mut exshift: f64 = 0.0;
    let mut p: f64 = 0.0;
    let mut q: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut z: f64 = 0.0;
    let mut t: f64 = 0.0;
    let mut w: f64 = 0.0;
    let mut x: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut norm: f64 = 0.0;
    for i in 0..nn {
        if i < low || i > high {
            d[i] = hmat[(i, i)];
            e[i] = 0.0;
        }
        for j in std::cmp::max(i - 1, 0)..nn {
            norm = norm + f64::abs(hmat[(i, j)]);
        }
    }
    let mut iter: isize = 0;
    while n >= low {
        let mut l: isize = n;
        while l > low {
            s = f64::abs(hmat[(l - 1, l - 1)]) + f64::abs(hmat[(l, l)]);
            if s == 0.0 {
                s = norm;
            }
            if f64::abs(hmat[(l, l - 1)]) <= eps * s {
                break;
            }
            l -= 1;
        }
        if l == n {
            hmat[(n, n)] = hmat[(n, n)] + exshift;
            d[n] = hmat[(n, n)];
            e[n] = 0.0;
            n -= 1;
            iter = 0;
        } else if l == n - 1 {
            w = hmat[(n, n - 1)] * hmat[(n - 1, n)];
            p = (hmat[(n - 1, n - 1)] - hmat[(n, n)]) / 2.0;
            q = p * p + w;
            z = f64::sqrt(f64::abs(q));
            hmat[(n, n)] = hmat[(n, n)] + exshift;
            hmat[(n - 1, n - 1)] = hmat[(n - 1, n - 1)] + exshift;
            x = hmat[(n, n)];
            if q >= 0.0 {
                if p >= 0.0 {
                    z = p + z;
                } else {
                    z = p - z;
                }
                d[n - 1] = x + z;
                d[n] = d[n - 1];
                if z != 0.0 {
                    d[n] = x - w / z;
                }
                e[n - 1] = 0.0;
                e[n] = 0.0;
                x = hmat[(n, n - 1)];
                s = f64::abs(x) + f64::abs(z);
                p = x / s;
                q = z / s;
                r = f64::sqrt(p * p + q * q);
                p = p / r;
                q = q / r;
                for j in n - 1..nn {
                    z = hmat[(n - 1, j)];
                    hmat[(n - 1, j)] = q * z + p * hmat[(n, j)];
                    hmat[(n, j)] = q * hmat[(n, j)] - p * z;
                }
                for i in 0..=n {
                    z = hmat[(i, n - 1)];
                    hmat[(i, n - 1)] = q * z + p * hmat[(i, n)];
                    hmat[(i, n)] = q * hmat[(i, n)] - p * z;
                }
                for i in low..=high {
                    z = vmat[(i, n - 1)];
                    vmat[(i, n - 1)] = q * z + p * vmat[(i, n)];
                    vmat[(i, n)] = q * vmat[(i, n)] - p * z;
                }
            } else {
                d[n - 1] = x + p;
                d[n] = x + p;
                e[n - 1] = z;
                e[n] = -z;
            }
            n = n - 2;
            iter = 0;
        } else {
            x = hmat[(n, n)];
            y = 0.0;
            w = 0.0;
            if l < n {
                y = hmat[(n - 1, n - 1)];
                w = hmat[(n, n - 1)] * hmat[(n - 1, n)];
            }
            if iter == 10 {
                exshift += x;
                for i in low..=n {
                    hmat[(i, i)] = hmat[(i, i)] - (x);
                }
                s = f64::abs(hmat[(n, n - 1)]) + f64::abs(hmat[(n - 1, n - 2)]);
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }
            if iter == 30 {
                s = (y - x) / 2.0;
                s = s * s + w;
                if s > 0.0 {
                    s = f64::sqrt(s);
                    if y < x {
                        s = -s;
                    }
                    s = x - w / ((y - x) / 2.0 + s);
                    for i in low..=n {
                        hmat[(i, i)] = hmat[(i, i)] - (s);
                    }
                    exshift += s;
                    x = 0.964;
                    y = x;
                    w = x;
                }
            }
            iter = iter + 1;
            if iter > max_iterations {
                return Err(MatrixError::new(
                    "general eigendecomposition did not converge within the iteration limit",
                ));
            }
            let mut m: isize = n - 2;
            while m >= l {
                z = hmat[(m, m)];
                r = x - z;
                s = y - z;
                p = (r * s - w) / hmat[(m + 1, m)] + hmat[(m, m + 1)];
                q = hmat[(m + 1, m + 1)] - z - r - s;
                r = hmat[(m + 2, m + 1)];
                s = f64::abs(p) + f64::abs(q) + f64::abs(r);
                p = p / s;
                q = q / s;
                r = r / s;
                if m == l {
                    break;
                }
                if f64::abs(hmat[(m, m - 1)]) * (f64::abs(q) + f64::abs(r))
                    < eps
                        * (f64::abs(p)
                            * (f64::abs(hmat[(m - 1, m - 1)])
                                + f64::abs(z)
                                + f64::abs(hmat[(m + 1, m + 1)])))
                {
                    break;
                }
                m -= 1;
            }
            for i in m + 2..=n {
                hmat[(i, i - 2)] = 0.0;
                if i > m + 2 {
                    hmat[(i, i - 3)] = 0.0;
                }
            }
            for k in m..=n - 1 {
                let notlast = k != n - 1;
                if k != m {
                    p = hmat[(k, k - 1)];
                    q = hmat[(k + 1, k - 1)];
                    r = 0.0;
                    if notlast {
                        r = hmat[(k + 2, k - 1)];
                    }
                    x = f64::abs(p) + f64::abs(q) + f64::abs(r);
                    if x == 0.0 {
                        continue;
                    }
                    p = p / x;
                    q = q / x;
                    r = r / x;
                }
                s = f64::sqrt(p * p + q * q + r * r);
                if p < 0.0 {
                    s = -s;
                }
                if s != 0.0 {
                    if k != m {
                        hmat[(k, k - 1)] = -s * x;
                    } else if l != m {
                        hmat[(k, k - 1)] = -hmat[(k, k - 1)];
                    }
                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;
                    for j in k..nn {
                        p = hmat[(k, j)] + q * hmat[(k + 1, j)];
                        if notlast {
                            p = p + r * hmat[(k + 2, j)];
                            hmat[(k + 2, j)] = hmat[(k + 2, j)] - p * z;
                        }
                        hmat[(k, j)] = hmat[(k, j)] - p * x;
                        hmat[(k + 1, j)] = hmat[(k + 1, j)] - p * y;
                    }
                    for i in 0..=std::cmp::min(n, k + 3) {
                        p = x * hmat[(i, k)] + y * hmat[(i, k + 1)];
                        if notlast {
                            p = p + z * hmat[(i, k + 2)];
                            hmat[(i, k + 2)] = hmat[(i, k + 2)] - p * r;
                        }
                        hmat[(i, k)] = hmat[(i, k)] - p;
                        hmat[(i, k + 1)] = hmat[(i, k + 1)] - p * q;
                    }
                    for i in low..=high {
                        p = x * vmat[(i, k)] + y * vmat[(i, k + 1)];
                        if notlast {
                            p = p + z * vmat[(i, k + 2)];
                            vmat[(i, k + 2)] = vmat[(i, k + 2)] - p * r;
                        }
                        vmat[(i, k)] = vmat[(i, k)] - p;
                        vmat[(i, k + 1)] = vmat[(i, k + 1)] - p * q;
                    }
                }
            }
        }
    }
    if norm == 0.0 {
        return Ok(());
    }
    for n in (0..=nn - 1).rev() {
        p = d[n];
        q = e[n];
        if q == 0.0 {
            let mut l: isize = n;
            hmat[(n, n)] = 1.0;
            for i in (0..=n - 1).rev() {
                w = hmat[(i, i)] - p;
                r = 0.0;
                for j in l..=n {
                    r = r + hmat[(i, j)] * hmat[(j, n)];
                }
                if e[i] < 0.0 {
                    z = w;
                    s = r;
                } else {
                    l = i;
                    if e[i] == 0.0 {
                        if w != 0.0 {
                            hmat[(i, n)] = -r / w;
                        } else {
                            hmat[(i, n)] = -r / (eps * norm);
                        }
                    } else {
                        x = hmat[(i, i + 1)];
                        y = hmat[(i + 1, i)];
                        q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                        t = (x * s - z * r) / q;
                        hmat[(i, n)] = t;
                        if f64::abs(x) > f64::abs(z) {
                            hmat[(i + 1, n)] = (-r - w * t) / x;
                        } else {
                            hmat[(i + 1, n)] = (-s - y * t) / z;
                        }
                    }
                    t = f64::abs(hmat[(i, n)]);
                    if (eps * t) * t > 1.0 {
                        for j in i..=n {
                            hmat[(j, n)] = hmat[(j, n)] / t;
                        }
                    }
                }
            }
        } else if q < 0.0 {
            let mut l: isize = n - 1;
            if f64::abs(hmat[(n, n - 1)]) > f64::abs(hmat[(n - 1, n)]) {
                hmat[(n - 1, n - 1)] = q / hmat[(n, n - 1)];
                hmat[(n - 1, n)] = -(hmat[(n, n)] - p) / hmat[(n, n - 1)];
            } else {
                (cdivr, cdivi) = general_cdiv(0.0, -hmat[(n - 1, n)], hmat[(n - 1, n - 1)] - p, q);
                hmat[(n - 1, n - 1)] = cdivr;
                hmat[(n - 1, n)] = cdivi;
            }
            hmat[(n, n - 1)] = 0.0;
            hmat[(n, n)] = 1.0;
            for i in (0..=n - 2).rev() {
                let mut ra: f64 = 0.0;
                let mut sa: f64 = 0.0;
                let mut vr: f64 = 0.0;
                let mut vi: f64 = 0.0;
                ra = 0.0;
                sa = 0.0;
                for j in l..=n {
                    ra = ra + hmat[(i, j)] * hmat[(j, n - 1)];
                    sa = sa + hmat[(i, j)] * hmat[(j, n)];
                }
                w = hmat[(i, i)] - p;
                if e[i] < 0.0 {
                    z = w;
                    r = ra;
                    s = sa;
                } else {
                    l = i;
                    if e[i] == 0.0 {
                        (cdivr, cdivi) = general_cdiv(-ra, -sa, w, q);
                        hmat[(i, n - 1)] = cdivr;
                        hmat[(i, n)] = cdivi;
                    } else {
                        x = hmat[(i, i + 1)];
                        y = hmat[(i + 1, i)];
                        vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                        vi = (d[i] - p) * 2.0 * q;
                        if vr == 0.0 && vi == 0.0 {
                            vr = eps
                                * norm
                                * (f64::abs(w)
                                    + f64::abs(q)
                                    + f64::abs(x)
                                    + f64::abs(y)
                                    + f64::abs(z));
                        }
                        (cdivr, cdivi) =
                            general_cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                        hmat[(i, n - 1)] = cdivr;
                        hmat[(i, n)] = cdivi;
                        if f64::abs(x) > (f64::abs(z) + f64::abs(q)) {
                            hmat[(i + 1, n - 1)] =
                                (-ra - w * hmat[(i, n - 1)] + q * hmat[(i, n)]) / x;
                            hmat[(i + 1, n)] = (-sa - w * hmat[(i, n)] - q * hmat[(i, n - 1)]) / x;
                        } else {
                            (cdivr, cdivi) = general_cdiv(
                                -r - y * hmat[(i, n - 1)],
                                -s - y * hmat[(i, n)],
                                z,
                                q,
                            );
                            hmat[(i + 1, n - 1)] = cdivr;
                            hmat[(i + 1, n)] = cdivi;
                        }
                    }
                    t = f64::max(f64::abs(hmat[(i, n - 1)]), f64::abs(hmat[(i, n)]));
                    if (eps * t) * t > 1.0 {
                        for j in i..=n {
                            hmat[(j, n - 1)] = hmat[(j, n - 1)] / t;
                            hmat[(j, n)] = hmat[(j, n)] / t;
                        }
                    }
                }
            }
        }
    }
    for i in 0..nn {
        if i < low || i > high {
            for j in i..nn {
                vmat[(i, j)] = hmat[(i, j)];
            }
        }
    }
    for j in (low..=nn - 1).rev() {
        for i in low..=high {
            z = 0.0;
            for k in low..=std::cmp::min(j, high) {
                z = z + vmat[(i, k)] * hmat[(k, j)];
            }
            vmat[(i, j)] = z;
        }
    }
    Ok(())
}

fn scale_binary(mut value: f64, mut power: i32) -> f64 {
    while power > 512 {
        value *= 2.0_f64.powi(512);
        power -= 512;
    }
    while power < -512 {
        value *= 2.0_f64.powi(-512);
        power += 512;
    }
    value * 2.0_f64.powi(power)
}
fn general_balance(h: &mut WorkMatrix) -> Vec<i32> {
    let n = h.n;
    let mut exponents = vec![0; n];
    fn log_norm(h: &WorkMatrix, i: usize, row: bool) -> f64 {
        let mut largest = 0.0_f64;
        for j in 0..h.n {
            if i != j {
                largest = largest.max(if row {
                    h.data[i * h.n + j].abs()
                } else {
                    h.data[j * h.n + i].abs()
                });
            }
        }
        if largest == 0.0 {
            return f64::NEG_INFINITY;
        }
        let mut sum = 0.0;
        for j in 0..h.n {
            if i != j {
                sum += if row {
                    h.data[i * h.n + j].abs() / largest
                } else {
                    h.data[j * h.n + i].abs() / largest
                };
            }
        }
        largest.log2() + sum.log2()
    }
    fn log_add(a: f64, b: f64) -> f64 {
        let top = a.max(b);
        top + ((a - top).exp2() + (b - top).exp2()).log2()
    }
    for _ in 0..64 {
        let mut changed = false;
        for i in 0..n {
            let r = log_norm(h, i, true);
            let c = log_norm(h, i, false);
            if !r.is_finite() || !c.is_finite() {
                continue;
            }
            let shift = (((r - c) / 2.0 + 0.5).floor() as i32).clamp(-512, 512);
            if shift == 0
                || log_add(r - shift as f64, c + shift as f64) >= log_add(r, c) + 0.95_f64.log2()
            {
                continue;
            }
            let mut safe = true;
            for j in 0..n {
                if i != j {
                    for (v, power) in [(h.data[i * n + j], -shift), (h.data[j * n + i], shift)] {
                        let scaled = scale_binary(v, power);
                        if !scaled.is_finite()
                            || (v != 0.0 && (scaled == 0.0 || scale_binary(scaled, -power) != v))
                        {
                            safe = false;
                        }
                    }
                }
            }
            if !safe {
                continue;
            }
            for j in 0..n {
                if i != j {
                    h.data[i * n + j] = scale_binary(h.data[i * n + j], -shift);
                    h.data[j * n + i] = scale_binary(h.data[j * n + i], shift);
                }
            }
            exponents[i] += shift;
            changed = true;
        }
        if !changed {
            break;
        }
    }
    exponents
}
fn binary_exponent(value: f64) -> i32 {
    let bits = (value.abs().to_bits() >> 52) & 0x7ff;
    if bits == 0 {
        ((value.abs() * 2.0_f64.powi(54)).to_bits() >> 52) as i32 - 1023 - 54
    } else {
        bits as i32 - 1023
    }
}
fn general_unbalance(vr: &mut Matrix, vi: &mut Matrix, col: usize, exponents: &[i32]) {
    let n = vr.rows;
    let mut common = i32::MIN;
    for row in 0..n {
        let v = vr.values[row * n + col]
            .abs()
            .max(vi.values[row * n + col].abs());
        if v != 0.0 {
            common = common.max(binary_exponent(v) + exponents[row]);
        }
    }
    if common == i32::MIN {
        return;
    }
    for row in 0..n {
        vr.values[row * n + col] = scale_binary(vr.values[row * n + col], exponents[row] - common);
        vi.values[row * n + col] = scale_binary(vi.values[row * n + col], exponents[row] - common);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn verify(a: &Matrix, g: &GeneralEigenDecomposition) {
        let n = a.rows;
        let mut scale = a.values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        if scale == 0.0 {
            scale = 1.0;
        }
        assert_eq!(g.eigenvalues_real.len(), n);
        assert_eq!(g.eigenvalues_imag.len(), n);
        assert_eq!((g.eigenvectors_real.rows, g.eigenvectors_imag.cols), (n, n));
        for j in 0..n {
            if j > 0 {
                assert!(
                    g.eigenvalues_real[j] > g.eigenvalues_real[j - 1]
                        || (g.eigenvalues_real[j] == g.eigenvalues_real[j - 1]
                            && g.eigenvalues_imag[j] >= g.eigenvalues_imag[j - 1])
                );
            }
            let (r, im) = (g.eigenvalues_real[j] / scale, g.eigenvalues_imag[j] / scale);
            let mut norm = 0.0;
            for i in 0..n {
                let (vr, vi) = (
                    g.eigenvectors_real.values[i * n + j],
                    g.eigenvectors_imag.values[i * n + j],
                );
                norm += vr * vr + vi * vi;
                let ar = (0..n)
                    .map(|k| a.values[i * n + k] / scale * g.eigenvectors_real.values[k * n + j])
                    .sum::<f64>();
                let ai = (0..n)
                    .map(|k| a.values[i * n + k] / scale * g.eigenvectors_imag.values[k * n + j])
                    .sum::<f64>();
                assert!(
                    (ar - r * vr + im * vi).hypot(ai - r * vi - im * vr) < 1e-8 * n.max(1) as f64,
                    "column {j}, row {i}"
                );
            }
            assert!((norm - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn complex_degenerate_and_extreme_scales() {
        for data in [
            vec![],
            vec![7.0],
            vec![0.0, -1.0, 1.0, 0.0],
            vec![2.0, 1.0, 0.0, 2.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0; 4],
            vec![1e300, 0.0, 0.0, 1e-300],
            vec![0.0, 1e300, -1e-300, 0.0],
        ] {
            let n = (data.len() as f64).sqrt() as usize;
            let a = Matrix::new(n, n, &data).unwrap();
            let g = a.eigen_general().unwrap();
            verify(&a, &g);
            assert_eq!(a.values, data);
            if data.len() == 4 && data[1] == 1e300 {
                assert!((g.eigenvalues_imag[0] + 1.0).abs() < 1e-12);
                assert!((g.eigenvalues_imag[1] - 1.0).abs() < 1e-12);
            }
        }
        for scale in [1e300, 1e-300, 1e-310] {
            let a = Matrix::new(2, 2, &[0.0, -scale, scale, 0.0]).unwrap();
            let g = a.eigen_general().unwrap();
            verify(&a, &g);
            assert!((g.eigenvalues_imag[1] / scale - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn dense_residuals_and_limits() {
        for n in 1..=12 {
            let values = (0..n * n)
                .map(|i| ((i * 17 + n * 11) % 31) as f64 / 8.0 - 15.0 / 8.0)
                .collect::<Vec<_>>();
            let a = Matrix::new(n, n, &values).unwrap();
            verify(&a, &a.eigen_general().unwrap());
        }
        let mut a = Matrix::new(
            4,
            4,
            &[
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 2.0, 5.0, 1.0, 3.0, 7.0, 2.0, 8.0, 1.0,
            ],
        )
        .unwrap();
        assert!(a.eigen_general_with(1).is_err());
        assert!(a.eigen_general_with(0).is_err());
        assert!(a.eigen_general_with(100001).is_err());
        assert!(Matrix::zeros(2, 3).unwrap().eigen_general().is_err());
        a.values[0] = f64::NAN;
        assert!(a.eigen_general().is_err());
    }
    #[test]
    fn triangular_extreme_diagonals_are_preserved() {
        for values in [[1e300, 1.0, 0.0, 1e-300], [1e300, 0.0, 1.0, 1e-300]] {
            let a = Matrix::new(2, 2, &values).unwrap();
            let g = a.eigen_general().unwrap();
            assert_eq!(g.eigenvalues_real, vec![1e-300, 1e300]);
            assert_eq!(g.eigenvalues_imag, vec![0.0, 0.0]);
            verify(&a, &g);
        }
        assert!(Matrix::new(2, 2, &[f64::MAX; 4])
            .unwrap()
            .eigen_general()
            .is_err());
    }
}
