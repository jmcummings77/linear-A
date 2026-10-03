use crate::{Matrix, MatrixError};

#[derive(Clone, Debug)]
pub struct CSRMatrix {
    rows: usize,
    cols: usize,
    rp: Vec<usize>,
    ci: Vec<usize>,
    values: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct CGResult {
    pub x: Vec<f64>,
    pub converged: bool,
    pub iterations: usize,
    pub reason: &'static str,
    pub residuals: Vec<f64>,
    pub iterates: Vec<Vec<f64>>,
}
#[derive(Clone, Copy, Debug)]
pub struct CGOptions {
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    pub max_iterations: usize,
    pub jacobi: bool,
    pub capture: bool,
}
impl Default for CGOptions {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            absolute_tolerance: 0.,
            max_iterations: 1000,
            jacobi: false,
            capture: false,
        }
    }
}
impl CSRMatrix {
    pub fn new(
        rows: usize,
        cols: usize,
        rp: &[usize],
        ci: &[usize],
        values: &[f64],
    ) -> Result<Self, MatrixError> {
        if rows.checked_add(1) != Some(rp.len())
            || ci.len() != values.len()
            || rp.first() != Some(&0)
            || rp.last() != Some(&values.len())
        {
            return Err(MatrixError::new("invalid CSR dimensions or arrays"));
        }
        if rp.iter().any(|&x| x > values.len()) || values.iter().any(|x| !x.is_finite()) {
            return Err(MatrixError::new("invalid CSR offsets or values"));
        }
        for i in 0..rows {
            if rp[i] > rp[i + 1] {
                return Err(MatrixError::new("CSR offsets must be monotone"));
            }
            for p in rp[i]..rp[i + 1] {
                if ci[p] >= cols || (p > rp[i] && ci[p] <= ci[p - 1]) {
                    return Err(MatrixError::new(
                        "CSR columns must be sorted unique and in range",
                    ));
                }
            }
        }
        Ok(Self {
            rows,
            cols,
            rp: rp.to_vec(),
            ci: ci.to_vec(),
            values: values.to_vec(),
        })
    }
    pub fn from_dense(a: &Matrix) -> Result<Self, MatrixError> {
        let mut rp = vec![0];
        let mut ci = vec![];
        let mut v = vec![];
        for i in 0..a.rows {
            for j in 0..a.cols {
                let x = a.values[i * a.cols + j];
                if x != 0. {
                    ci.push(j);
                    v.push(x);
                }
            }
            rp.push(v.len());
        }
        Self::new(a.rows, a.cols, &rp, &ci, &v)
    }
    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
    pub fn row_offsets(&self) -> &[usize] {
        &self.rp
    }
    pub fn column_indices(&self) -> &[usize] {
        &self.ci
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn matvec(&self, x: &[f64]) -> Result<Vec<f64>, MatrixError> {
        if x.len() != self.cols || x.iter().any(|v| !v.is_finite()) {
            return Err(MatrixError::new("invalid vector"));
        }
        (0..self.rows)
            .map(|i| {
                let sum = (self.rp[i]..self.rp[i + 1])
                    .map(|p| self.values[p] * x[self.ci[p]])
                    .sum::<f64>();
                crate::finite(sum)
            })
            .collect()
    }
    fn find(&self, row: usize, col: usize) -> usize {
        let mut lo = self.rp[row];
        let mut hi = self.rp[row + 1];
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.ci[mid] < col {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }
    pub fn conjugate_gradient(&self, b: &[f64], o: CGOptions) -> Result<CGResult, MatrixError> {
        let n = self.rows;
        if self.cols != n || b.len() != n || b.iter().any(|x| !x.is_finite()) {
            return Err(MatrixError::new(
                "CG requires square matrix and finite matching vector",
            ));
        }
        if !o.relative_tolerance.is_finite()
            || o.relative_tolerance < 0.
            || o.relative_tolerance >= 1.
            || !o.absolute_tolerance.is_finite()
            || o.absolute_tolerance < 0.
            || o.max_iterations > 100000
        {
            return Err(MatrixError::new("invalid CG options"));
        }
        let mut diagonal = vec![1.; n];
        for i in 0..n {
            for p in self.rp[i]..self.rp[i + 1] {
                let j = self.ci[p];
                let q = self.find(j, i);
                let other = if q < self.rp[j + 1] && self.ci[q] == i {
                    self.values[q]
                } else {
                    0.
                };
                if self.values[p] != other {
                    return Err(MatrixError::new("CG requires exact symmetry"));
                }
            }
            if o.jacobi {
                let q = self.find(i, i);
                if q == self.rp[i + 1] || self.ci[q] != i || self.values[q] <= 0. {
                    return Err(MatrixError::new("Jacobi requires positive diagonal"));
                }
                diagonal[i] = self.values[q];
            }
        }
        let norm = |v: &[f64]| v.iter().fold(0_f64, |s, &x| s.hypot(x));
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let mut x = vec![0.; n];
        let mut r = b.to_vec();
        let mut history = vec![norm(&r)];
        let mut frames = if o.capture { vec![x.clone()] } else { vec![] };
        let threshold = o.absolute_tolerance.max(o.relative_tolerance * history[0]);
        let mut reason = if !history[0].is_finite() {
            "nonfinite"
        } else if history[0] <= threshold {
            "converged"
        } else {
            "iteration_limit"
        };
        let mut z: Vec<f64> = r.iter().enumerate().map(|(i, v)| v / diagonal[i]).collect();
        let mut p = z.clone();
        let mut rho = dot(&r, &z);
        if reason == "iteration_limit" {
            for _ in 0..o.max_iterations {
                let q = match self.matvec(&p) {
                    Ok(v) => v,
                    Err(_) => {
                        reason = "nonfinite";
                        break;
                    }
                };
                let curvature = dot(&p, &q);
                if !rho.is_finite() || !curvature.is_finite() {
                    reason = "nonfinite";
                    break;
                }
                if rho <= 0. || curvature <= 0. {
                    reason = "breakdown";
                    break;
                }
                let alpha = rho / curvature;
                let candidate: Vec<f64> = x
                    .iter()
                    .enumerate()
                    .map(|(i, v)| v + alpha * p[i])
                    .collect();
                let ax = match self.matvec(&candidate) {
                    Ok(v) => v,
                    Err(_) => {
                        reason = "nonfinite";
                        break;
                    }
                };
                let residual: Vec<f64> = b.iter().zip(ax).map(|(v, y)| v - y).collect();
                let length = norm(&residual);
                if !length.is_finite() {
                    reason = "nonfinite";
                    break;
                }
                x = candidate;
                r = residual;
                history.push(length);
                if o.capture {
                    frames.push(x.clone());
                }
                if length <= threshold {
                    reason = "converged";
                    break;
                }
                z = r.iter().enumerate().map(|(i, v)| v / diagonal[i]).collect();
                let next = dot(&r, &z);
                if !next.is_finite() {
                    reason = "nonfinite";
                    break;
                }
                if next <= 0. {
                    reason = "breakdown";
                    break;
                }
                let beta = next / rho;
                p = z.iter().enumerate().map(|(i, v)| v + beta * p[i]).collect();
                rho = next;
            }
        }
        Ok(CGResult {
            x,
            converged: reason == "converged",
            iterations: history.len() - 1,
            reason,
            residuals: history,
            iterates: frames,
        })
    }
}
