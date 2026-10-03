use crate::{Matrix, MatrixError};

/// Economy factors A=U diag(values) Vt, with descending nonnegative values.
#[derive(Clone, Debug)]
pub struct SingularValueDecomposition {
    pub u: Matrix,
    pub values: Vec<f64>,
    pub vt: Matrix,
}
impl Matrix {
    pub fn svd(&self) -> Result<SingularValueDecomposition, MatrixError> {
        self.svd_with(1e-12, 100)
    }
    /// Cyclic one-sided Jacobi; see ports/SVD.md for range and convergence rules.
    pub fn svd_with(
        &self,
        tolerance: f64,
        max_sweeps: usize,
    ) -> Result<SingularValueDecomposition, MatrixError> {
        let fail = |s| MatrixError::new(s);
        if !tolerance.is_finite()
            || tolerance <= 0.
            || tolerance >= 1.
            || max_sweeps == 0
            || max_sweeps > 10000
        {
            return Err(fail("invalid SVD options"));
        }
        let (m, n) = (self.rows, self.cols);
        if m < n {
            let r = self.transpose()?.svd_with(tolerance, max_sweeps)?;
            return Ok(SingularValueDecomposition {
                u: r.vt.transpose()?,
                values: r.values,
                vt: r.u.transpose()?,
            });
        }
        let mut scale: f64 = 0.;
        for &x in &self.values {
            if !x.is_finite() {
                return Err(fail("SVD requires finite input"));
            }
            scale = scale.max(x.abs());
        }
        if scale == 0. {
            scale = 1.;
        }
        let mut b = self.values.clone();
        for (i, x) in b.iter_mut().enumerate() {
            *x /= scale;
            if self.values[i] != 0. && *x == 0. {
                return Err(fail("SVD scaling discards an entry"));
            }
        }
        let mut v = Matrix::identity(n)?.values;
        let norm = |b: &[f64], j: usize| {
            let mut s: f64 = 0.;
            for i in 0..m {
                s = s.hypot(b[i * n + j]);
            }
            s
        };
        let mut converged = false;
        for sweep in 0..=max_sweeps {
            let mut changed = false;
            for p in 0..n {
                for q in p + 1..n {
                    let (np, nq) = (norm(&b, p), norm(&b, q));
                    if np == 0. || nq == 0. {
                        continue;
                    }
                    let mut corr = 0.;
                    for i in 0..m {
                        corr += (b[i * n + p] / np) * (b[i * n + q] / nq);
                    }
                    if corr.abs() <= tolerance {
                        continue;
                    }
                    changed = true;
                    if sweep == max_sweeps {
                        continue;
                    }
                    let pair = np.max(nq);
                    let (ap, aq) = (np / pair, nq / pair);
                    let (delta, g) = (aq * aq - ap * ap, 2. * ap * aq * corr);
                    let t = if delta == 0. {
                        1_f64.copysign(g)
                    } else {
                        g / (delta + delta.hypot(g).copysign(delta))
                    };
                    if t.abs() < f64::MIN_POSITIVE {
                        let small = if np < nq { p } else { q };
                        for i in 0..m {
                            b[i * n + small] = 0.;
                        }
                        continue;
                    }
                    let c = 1. / 1_f64.hypot(t);
                    let s = c * t;
                    for i in 0..m {
                        let (x, y) = (b[i * n + p], b[i * n + q]);
                        b[i * n + p] = c * x - s * y;
                        b[i * n + q] = s * x + c * y;
                    }
                    for i in 0..n {
                        let (x, y) = (v[i * n + p], v[i * n + q]);
                        v[i * n + p] = c * x - s * y;
                        v[i * n + q] = s * x + c * y;
                    }
                }
            }
            if !changed {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(fail("SVD did not converge"));
        }
        let norms: Vec<f64> = (0..n).map(|j| norm(&b, j)).collect();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| norms[b].partial_cmp(&norms[a]).unwrap());
        let mut u = Matrix::zeros(m, n)?;
        let mut vt = Matrix::zeros(n, n)?;
        let mut values = Vec::with_capacity(n);
        for (j, &k) in order.iter().enumerate() {
            let value = norms[k] * scale;
            if !value.is_finite() || (norms[k] != 0. && value == 0.) {
                return Err(fail("singular value outside float64 range"));
            }
            values.push(value);
            for i in 0..n {
                vt.values[j * n + i] = v[i * n + k];
            }
            if norms[k] != 0. {
                for i in 0..m {
                    u.values[i * n + j] = b[i * n + k] / norms[k];
                }
            } else {
                let mut found = false;
                for axis in 0..m {
                    let mut candidate = vec![0.; m];
                    candidate[axis] = 1.;
                    for _ in 0..2 {
                        for col in 0..j {
                            let mut dot = 0.;
                            for i in 0..m {
                                dot += candidate[i] * u.values[i * n + col];
                            }
                            for i in 0..m {
                                candidate[i] -= dot * u.values[i * n + col];
                            }
                        }
                    }
                    let mut length: f64 = 0.;
                    for &x in &candidate {
                        length = length.hypot(x);
                    }
                    if length > 0.5 / (m as f64).sqrt() {
                        for i in 0..m {
                            u.values[i * n + j] = candidate[i] / length;
                        }
                        found = true;
                        break;
                    }
                }
                if !found {
                    return Err(fail("cannot complete SVD null basis"));
                }
            }
        }
        Ok(SingularValueDecomposition { u, values, vt })
    }
}

/// Rank and reciprocal 2-norm conditions before and after truncation.
#[derive(Clone, Debug)]
pub struct SpectralDiagnostics {
    pub rank: usize,
    pub reciprocal_condition: f64,
    pub retained_reciprocal_condition: f64,
}
fn spectral_rank(s: &[f64], cutoff: f64) -> usize {
    s.iter()
        .take_while(|&&x| x > 0. && (cutoff == 0. || x / s[0] > cutoff))
        .count()
}
fn check_cutoff(cutoff: f64) -> Result<(), MatrixError> {
    if !cutoff.is_finite() || !(0. ..=1.).contains(&cutoff) {
        Err(MatrixError::new("invalid relative cutoff"))
    } else {
        Ok(())
    }
}
impl Matrix {
    pub fn spectral_diagnostics(&self) -> Result<SpectralDiagnostics, MatrixError> {
        self.spectral_diagnostics_with(self.rows.max(self.cols) as f64 * f64::EPSILON)
    }
    pub fn spectral_diagnostics_with(
        &self,
        cutoff: f64,
    ) -> Result<SpectralDiagnostics, MatrixError> {
        check_cutoff(cutoff)?;
        let r = self.svd()?;
        let s = &r.values;
        let rank = spectral_rank(s, cutoff);
        Ok(SpectralDiagnostics {
            rank,
            reciprocal_condition: if !s.is_empty() && s[0] > 0. {
                s[s.len() - 1] / s[0]
            } else {
                0.
            },
            retained_reciprocal_condition: if rank > 0 { s[rank - 1] / s[0] } else { 0. },
        })
    }
    pub fn pseudoinverse(&self) -> Result<Matrix, MatrixError> {
        self.pseudoinverse_with(self.rows.max(self.cols) as f64 * f64::EPSILON)
    }
    pub fn pseudoinverse_with(&self, cutoff: f64) -> Result<Matrix, MatrixError> {
        self.apply_inverse(None, cutoff, 0.)
    }
    pub fn solve_minimum_norm(&self, rhs: &Matrix) -> Result<Matrix, MatrixError> {
        self.solve_minimum_norm_with(rhs, self.rows.max(self.cols) as f64 * f64::EPSILON)
    }
    pub fn solve_minimum_norm_with(
        &self,
        rhs: &Matrix,
        cutoff: f64,
    ) -> Result<Matrix, MatrixError> {
        self.apply_inverse(Some(rhs), cutoff, 0.)
    }
    pub fn solve_ridge(&self, rhs: &Matrix, lambda: f64) -> Result<Matrix, MatrixError> {
        if !lambda.is_finite() || lambda < 0. {
            return Err(MatrixError::new("lambda must be finite and nonnegative"));
        }
        let cutoff = if lambda == 0. {
            self.rows.max(self.cols) as f64 * f64::EPSILON
        } else {
            0.
        };
        self.apply_inverse(Some(rhs), cutoff, lambda)
    }
    fn apply_inverse(
        &self,
        rhs: Option<&Matrix>,
        cutoff: f64,
        lambda: f64,
    ) -> Result<Matrix, MatrixError> {
        use crate::determinant::quotient_product;
        check_cutoff(cutoff)?;
        if rhs.map_or(false, |b| b.rows != self.rows) {
            return Err(MatrixError::new("incompatible right-hand side"));
        }
        let r = self.svd()?;
        let (m, n, k) = (self.rows, self.cols, r.values.len());
        let rank = spectral_rank(&r.values, cutoff);
        let cols = rhs.map_or(m, |b| b.cols);
        let mut out = Matrix::zeros(n, cols)?;
        for j in 0..cols {
            let scale = rhs.map_or(1., |b| {
                (0..m).fold(0_f64, |s, i| s.max(b.values[i * cols + j].abs()))
            });
            if rank > 0
                && scale != 0.
                && rhs.map_or(false, |b| {
                    (0..m).any(|i| {
                        b.values[i * cols + j] != 0. && b.values[i * cols + j] / scale == 0.
                    })
                })
            {
                return Err(MatrixError::new(
                    "right-hand side scaling discards an entry",
                ));
            }
            for p in 0..rank {
                let projection =
                    rhs.map_or(r.u.values.get(j * k + p).copied().unwrap_or(0.), |b| {
                        if scale == 0. {
                            0.
                        } else {
                            (0..m)
                                .map(|i| r.u.values[i * k + p] * (b.values[i * cols + j] / scale))
                                .sum()
                        }
                    });
                let coefficient = if rhs.is_some() {
                    crate::finite(if lambda > 0. {
                        crate::determinant::ridge_product(projection, r.values[p], scale, lambda)
                    } else {
                        quotient_product(projection, r.values[p], scale)
                    })?
                } else {
                    0.
                };
                for i in 0..n {
                    let term = if rhs.is_some() {
                        r.vt.values[p * n + i] * coefficient
                    } else {
                        quotient_product(projection, r.values[p], r.vt.values[p * n + i])
                    };
                    out.values[i * cols + j] = crate::finite(out.values[i * cols + j] + term)?;
                }
            }
        }
        Ok(out)
    }
}
