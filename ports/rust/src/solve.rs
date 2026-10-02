use crate::{allocate, finite, Matrix, MatrixError};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FactorAlgorithm {
    Lu,
    Cholesky,
    Qr,
}

/// Independent reusable factors. QR requires rows >= columns and numerical full rank.
#[derive(Clone, Debug)]
pub struct Factorization {
    rows: usize,
    cols: usize,
    algorithm: FactorAlgorithm,
    scale: f64,
    norm: f64,
    data: Vec<f64>,
    tau: Vec<f64>,
    permutation: Vec<usize>,
}
fn scaled(value: f64, scale: f64) -> Result<f64, MatrixError> {
    let result = finite(value / scale)?;
    if value != 0.0 && result == 0.0 {
        return Err(MatrixError::new(
            "solver scaling would discard a nonzero value",
        ));
    }
    Ok(result)
}
impl Matrix {
    pub fn factor_lu(&self) -> Result<Factorization, MatrixError> {
        Factorization::new(self, FactorAlgorithm::Lu)
    }
    pub fn factor_cholesky(&self) -> Result<Factorization, MatrixError> {
        Factorization::new(self, FactorAlgorithm::Cholesky)
    }
    pub fn factor_qr(&self) -> Result<Factorization, MatrixError> {
        Factorization::new(self, FactorAlgorithm::Qr)
    }
    pub fn solve(&self, rhs: &Matrix) -> Result<Matrix, MatrixError> {
        self.factor_lu()?.solve(rhs)
    }
    pub fn least_squares(&self, rhs: &Matrix) -> Result<Matrix, MatrixError> {
        self.factor_qr()?.solve(rhs)
    }
}
impl Factorization {
    pub fn new(source: &Matrix, algorithm: FactorAlgorithm) -> Result<Self, MatrixError> {
        let (m, n) = (source.rows, source.cols);
        if m < n || (algorithm != FactorAlgorithm::Qr && m != n) {
            return Err(MatrixError::new(
                "factorization requires square input, or rows >= columns for QR",
            ));
        }
        let maximum = source.values.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        // Extract a power of two directly, including subnormal inputs.
        let scale = if maximum == 0.0 {
            1.0
        } else {
            let bits = maximum.to_bits();
            let exponent = bits & (0x7ff_u64 << 52);
            if exponent != 0 {
                f64::from_bits(exponent)
            } else {
                f64::from_bits(1_u64 << (63 - bits.leading_zeros()))
            }
        };
        let mut a = allocate(source.values.len())?;
        for (out, &value) in a.iter_mut().zip(&source.values) {
            *out = scaled(value, scale)?;
        }
        let mut norm = 0.0_f64;
        for i in 0..m {
            norm = norm.max((0..n).map(|j| a[i * n + j].abs()).sum());
        }
        let mut permutation = Vec::new();
        permutation
            .try_reserve_exact(m)
            .map_err(|_| MatrixError::new("unable to allocate factor permutation"))?;
        permutation.extend(0..m);
        let mut tau = allocate(n)?;
        match algorithm {
            FactorAlgorithm::Lu => {
                for k in 0..n {
                    let mut pivot = k;
                    for i in k + 1..n {
                        if a[i * n + k].abs() > a[pivot * n + k].abs() {
                            pivot = i;
                        }
                    }
                    if a[pivot * n + k] == 0.0 {
                        return Err(MatrixError::new("singular matrix: zero computed LU pivot"));
                    }
                    for j in 0..n {
                        a.swap(k * n + j, pivot * n + j);
                    }
                    permutation.swap(k, pivot);
                    for i in k + 1..n {
                        a[i * n + k] = finite(a[i * n + k] / a[k * n + k])?;
                        for j in k + 1..n {
                            a[i * n + j] = finite(a[i * n + j] - a[i * n + k] * a[k * n + j])?;
                        }
                    }
                }
            }
            FactorAlgorithm::Cholesky => {
                for i in 0..n {
                    for j in 0..i {
                        if source.values[i * n + j] != source.values[j * n + i] {
                            return Err(MatrixError::new("Cholesky requires exact symmetry"));
                        }
                    }
                }
                for i in 0..n {
                    for j in 0..=i {
                        let mut value = a[i * n + j];
                        for k in 0..j {
                            value = finite(value - a[i * n + k] * a[j * n + k])?;
                        }
                        a[i * n + j] = if i == j {
                            if value <= 0.0 {
                                return Err(MatrixError::new(
                                    "Cholesky requires positive computed pivots",
                                ));
                            }
                            value.sqrt()
                        } else {
                            finite(value / a[j * n + j])?
                        };
                    }
                }
            }
            FactorAlgorithm::Qr => {
                let mut largest = 0.0_f64;
                for j in 0..n {
                    let mut norm = 0.0_f64;
                    for i in 0..m {
                        norm = norm.hypot(a[i * n + j]);
                    }
                    largest = largest.max(norm);
                }
                let threshold = f64::EPSILON * m as f64 * largest;
                for k in 0..n {
                    let (mut pivot, mut norm) = (k, 0.0_f64);
                    for j in k..n {
                        let mut candidate = 0.0_f64;
                        for i in k..m {
                            candidate = candidate.hypot(a[i * n + j]);
                        }
                        if candidate > norm {
                            pivot = j;
                            norm = candidate;
                        }
                    }
                    if norm <= threshold {
                        return Err(MatrixError::new("QR input is numerically rank deficient"));
                    }
                    for i in 0..m {
                        a.swap(i * n + k, i * n + pivot);
                    }
                    permutation.swap(k, pivot);
                    let old = a[k * n + k];
                    let alpha = -norm.copysign(old);
                    let divisor = old - alpha;
                    tau[k] = (alpha - old) / alpha;
                    for i in k + 1..m {
                        a[i * n + k] /= divisor;
                    }
                    a[k * n + k] = alpha;
                    for j in k + 1..n {
                        let mut dot = a[k * n + j];
                        for i in k + 1..m {
                            dot += a[i * n + k] * a[i * n + j];
                        }
                        dot *= tau[k];
                        a[k * n + j] = finite(a[k * n + j] - dot)?;
                        for i in k + 1..m {
                            a[i * n + j] = finite(a[i * n + j] - a[i * n + k] * dot)?;
                        }
                    }
                }
            }
        }
        Ok(Self {
            rows: m,
            cols: n,
            algorithm,
            scale,
            norm,
            data: a,
            tau,
            permutation,
        })
    }
    fn solve_internal(&self, rhs: &Matrix, rescale: bool) -> Result<Matrix, MatrixError> {
        let (m, n, p) = (self.rows, self.cols, rhs.cols);
        if rhs.rows != m {
            return Err(MatrixError::new(
                "right-hand side row count must match factorization",
            ));
        }
        let mut work = allocate(rhs.values.len())?;
        let a = &self.data;
        for i in 0..m {
            let row = if self.algorithm == FactorAlgorithm::Lu {
                self.permutation[i]
            } else {
                i
            };
            for j in 0..p {
                let value = rhs.values[row * p + j];
                work[i * p + j] = if rescale {
                    scaled(value, self.scale)?
                } else {
                    value
                };
            }
        }
        if self.algorithm == FactorAlgorithm::Qr {
            for k in 0..n {
                for j in 0..p {
                    let mut dot = work[k * p + j];
                    for i in k + 1..m {
                        dot = finite(dot + a[i * n + k] * work[i * p + j])?;
                    }
                    dot = finite(dot * self.tau[k])?;
                    work[k * p + j] = finite(work[k * p + j] - dot)?;
                    for i in k + 1..m {
                        work[i * p + j] = finite(work[i * p + j] - a[i * n + k] * dot)?;
                    }
                }
            }
        } else {
            for i in 0..n {
                for j in 0..p {
                    let mut value = work[i * p + j];
                    for k in 0..i {
                        value = finite(value - a[i * n + k] * work[k * p + j])?;
                    }
                    work[i * p + j] = if self.algorithm == FactorAlgorithm::Cholesky {
                        finite(value / a[i * n + i])?
                    } else {
                        value
                    };
                }
            }
        }
        for i in (0..n).rev() {
            for j in 0..p {
                let mut value = work[i * p + j];
                for k in i + 1..n {
                    let coefficient = if self.algorithm == FactorAlgorithm::Cholesky {
                        a[k * n + i]
                    } else {
                        a[i * n + k]
                    };
                    value = finite(value - coefficient * work[k * p + j])?;
                }
                work[i * p + j] = finite(value / a[i * n + i])?;
            }
        }
        let mut result = Matrix::zeros(n, p)?;
        for i in 0..n {
            let row = if self.algorithm == FactorAlgorithm::Qr {
                self.permutation[i]
            } else {
                i
            };
            result.values[row * p..(row + 1) * p].copy_from_slice(&work[i * p..(i + 1) * p]);
        }
        Ok(result)
    }
    pub fn solve(&self, rhs: &Matrix) -> Result<Matrix, MatrixError> {
        self.solve_internal(rhs, true)
    }
    /// Computed infinity-norm rcond using an inverse: O(n^3), not a certified bound.
    pub fn reciprocal_condition(&self) -> Result<f64, MatrixError> {
        let n = self.cols;
        if self.rows != n {
            return Err(MatrixError::new(
                "condition diagnostic requires a square matrix",
            ));
        }
        if n == 0 {
            return Ok(1.0);
        }
        let identity = Matrix::identity(n)?;
        let inverse = match self.solve_internal(&identity, false) {
            Ok(value) => value,
            Err(error) if error.0 == "matrix arithmetic produced a nonfinite value" => {
                return Ok(0.0)
            }
            Err(error) => return Err(error),
        };
        let mut norm = 0.0_f64;
        for i in 0..n {
            norm = norm.max((0..n).map(|j| inverse.values[i * n + j].abs()).sum());
        }
        Ok(((1.0 / self.norm) / norm).min(1.0))
    }
}
