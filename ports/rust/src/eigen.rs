use crate::{finite, Matrix, MatrixError};

/// Ascending real eigenvalues and the corresponding orthonormal eigenvector columns.
/// Repeated eigenvalues may have any orthonormal basis of their eigenspace.
#[derive(Clone, Debug)]
pub struct SymmetricEigenDecomposition {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Matrix,
}

impl Matrix {
    /// Solve a finite, exactly symmetric matrix using cyclic Jacobi rotations.
    /// Uses a relative Frobenius tolerance of 1e-12 and at most 50 sweeps.
    /// Does not modify the input. Tiny eigenvalues in mixed-scale matrices may
    /// have large relative error even when the normwise residual is small.
    pub fn eigen_symmetric(&self) -> Result<SymmetricEigenDecomposition, MatrixError> {
        self.eigen_symmetric_with(1e-12, 50)
    }

    /// Set the relative tolerance (finite and between zero and one) and positive
    /// sweep limit. Returns an error on nonconvergence or nonfinite eigenvalues.
    pub fn eigen_symmetric_with(
        &self,
        tolerance: f64,
        max_sweeps: usize,
    ) -> Result<SymmetricEigenDecomposition, MatrixError> {
        self.require_square()?;
        if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0 || max_sweeps == 0 {
            return Err(MatrixError::new(
                "invalid eigensolver tolerance or sweep limit",
            ));
        }
        let n = self.rows;
        let mut scale: f64 = 0.0;
        let mut diagonal = true;
        for row in 0..n {
            for col in 0..n {
                let value = self.values[row * n + col];
                if !value.is_finite() || value != self.values[col * n + row] {
                    return Err(MatrixError::new(
                        "eigendecomposition requires a finite exactly symmetric matrix",
                    ));
                }
                scale = scale.max(value.abs());
                if row != col && value != 0.0 {
                    diagonal = false;
                }
            }
        }
        let mut vectors = Matrix::identity(n)?;
        let mut work = self.values.clone();
        if diagonal {
            return finish(&work, vectors, 1.0);
        }
        let exponent_bits = (scale.to_bits() >> 52) & 0x7ff;
        if exponent_bits != 0 {
            scale = f64::from_bits(exponent_bits << 52);
        }
        let mut norm: f64 = 0.0;
        for value in &mut work {
            *value /= scale;
            norm = norm.hypot(*value);
        }
        let threshold = tolerance * norm;
        let pair_threshold = threshold / (2.0 * n as f64);
        for sweep in 0..=max_sweeps {
            let mut off_norm: f64 = 0.0;
            for p in 0..n {
                for q in p + 1..n {
                    off_norm = off_norm.hypot(std::f64::consts::SQRT_2 * work[p * n + q]);
                }
            }
            if off_norm <= threshold {
                return finish(&work, vectors, scale);
            }
            if sweep == max_sweeps {
                break;
            }
            for p in 0..n {
                for q in p + 1..n {
                    let b = work[p * n + q];
                    if b.abs() <= pair_threshold {
                        continue;
                    }
                    let delta = (work[q * n + q] - work[p * n + p]) / 2.0;
                    let hypotenuse = delta.hypot(b);
                    let t = b / (delta + if delta < 0.0 { -hypotenuse } else { hypotenuse });
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    work[p * n + p] -= t * b;
                    work[q * n + q] += t * b;
                    work[p * n + q] = 0.0;
                    work[q * n + p] = 0.0;
                    for k in 0..n {
                        if k != p && k != q {
                            let x = work[k * n + p];
                            let y = work[k * n + q];
                            let xp = c * x - s * y;
                            let yq = s * x + c * y;
                            work[k * n + p] = xp;
                            work[p * n + k] = xp;
                            work[k * n + q] = yq;
                            work[q * n + k] = yq;
                        }
                        let vp = vectors.values[k * n + p];
                        let vq = vectors.values[k * n + q];
                        vectors.values[k * n + p] = c * vp - s * vq;
                        vectors.values[k * n + q] = s * vp + c * vq;
                    }
                }
            }
        }
        Err(MatrixError::new(
            "symmetric eigendecomposition did not converge within the sweep limit",
        ))
    }
}

fn finish(
    work: &[f64],
    vectors: Matrix,
    scale: f64,
) -> Result<SymmetricEigenDecomposition, MatrixError> {
    let n = vectors.rows;
    let values: Vec<f64> = (0..n)
        .map(|i| finite(work[i * n + i] * scale))
        .collect::<Result<_, _>>()?;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| values[a].partial_cmp(&values[b]).unwrap().then(a.cmp(&b)));
    let mut sorted = Matrix::zeros(n, n)?;
    for (col, &original) in order.iter().enumerate() {
        let mut length: f64 = 0.0;
        let mut largest = 0;
        for row in 0..n {
            let value = vectors.values[row * n + original];
            length = length.hypot(value);
            if value.abs() > vectors.values[largest * n + original].abs() {
                largest = row;
            }
        }
        if vectors.values[largest * n + original] < 0.0 {
            length = -length;
        }
        for row in 0..n {
            sorted.values[row * n + col] = vectors.values[row * n + original] / length;
        }
    }
    Ok(SymmetricEigenDecomposition {
        eigenvalues: order.iter().map(|&i| values[i]).collect(),
        eigenvectors: sorted,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn verify(a: &Matrix, result: &SymmetricEigenDecomposition) {
        let n = a.rows;
        let scale = a
            .values
            .iter()
            .fold(0.0_f64, |acc, x| acc.max(x.abs()))
            .max(f64::MIN_POSITIVE);
        let q = &result.eigenvectors.values;
        for i in 0..n {
            for j in 0..n {
                let aq: f64 = (0..n)
                    .map(|k| (a.values[i * n + k] / scale) * q[k * n + j])
                    .sum();
                let dot: f64 = (0..n).map(|k| q[k * n + i] * q[k * n + j]).sum();
                assert!((aq - q[i * n + j] * (result.eigenvalues[j] / scale)).abs() < 1e-9);
                assert!((dot - if i == j { 1.0 } else { 0.0 }).abs() < 1e-10);
            }
        }
        assert!(result.eigenvalues.windows(2).all(|p| p[0] <= p[1]));
    }

    #[test]
    fn analytic_and_dense_residuals() {
        for n in 0..=12 {
            let mut a = Matrix::zeros(n, n).unwrap();
            for i in 0..n {
                for j in i..n {
                    let value = ((i * 17 + j * 13) % 19) as f64 / 8.0 - 1.0;
                    a.values[i * n + j] = value;
                    a.values[j * n + i] = value;
                }
            }
            let original = a.clone();
            verify(&a, &a.eigen_symmetric().unwrap());
            assert_eq!(a, original);
        }
        for scale in [1.0, 1e300, 1e-300, 1e-310] {
            let a = Matrix::new(2, 2, &[2.0 * scale, scale, scale, 2.0 * scale]).unwrap();
            let result = a.eigen_symmetric().unwrap();
            assert!((result.eigenvalues[0] / scale - 1.0).abs() < 1e-12);
            assert!((result.eigenvalues[1] / scale - 3.0).abs() < 1e-12);
            verify(&a, &result);
        }
    }

    #[test]
    fn diagonal_repeated_and_tight_tolerance() {
        let diagonal =
            Matrix::new(3, 3, &[1e300, 0.0, 0.0, 0.0, 1e-300, 0.0, 0.0, 0.0, -2.0]).unwrap();
        assert_eq!(
            diagonal.eigen_symmetric().unwrap().eigenvalues,
            vec![-2.0, 1e-300, 1e300]
        );
        verify(
            &Matrix::identity(3).unwrap(),
            &Matrix::identity(3).unwrap().eigen_symmetric().unwrap(),
        );
        let tiny = Matrix::new(2, 2, &[1.0, 1e-200, 1e-200, 2.0])
            .unwrap()
            .eigen_symmetric_with(1e-250, 1)
            .unwrap();
        assert!((tiny.eigenvectors.values[1] / 1e-200 - 1.0).abs() < 1e-14);
    }

    #[test]
    fn validation_and_nonconvergence() {
        assert!(Matrix::zeros(2, 3).unwrap().eigen_symmetric().is_err());
        assert!(Matrix::new(2, 2, &[1.0, 1.0, 0.0, 1.0])
            .unwrap()
            .eigen_symmetric()
            .is_err());
        let a = Matrix::new(3, 3, &[1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0]).unwrap();
        for tolerance in [0.0, -1.0, 1.0, f64::NAN, f64::INFINITY] {
            assert!(a.eigen_symmetric_with(tolerance, 50).is_err());
        }
        assert!(a.eigen_symmetric_with(1e-12, 0).is_err());
        assert!(a.eigen_symmetric_with(1e-12, 1).is_err());
        assert!(Matrix::new(2, 2, &[f64::MAX; 4])
            .unwrap()
            .eigen_symmetric()
            .is_err());
    }
}
