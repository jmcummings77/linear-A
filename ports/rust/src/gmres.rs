use crate::{CSRMatrix, MatrixError};
#[derive(Clone, Copy, Debug)]
pub struct GMRESOptions {
    pub restart: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    pub max_iterations: usize,
    pub jacobi: bool,
    pub capture: bool,
}
impl Default for GMRESOptions {
    fn default() -> Self {
        Self {
            restart: 30,
            relative_tolerance: 1e-10,
            absolute_tolerance: 0.,
            max_iterations: 1000,
            jacobi: false,
            capture: false,
        }
    }
}
#[derive(Clone, Debug)]
pub struct GMRESResult {
    pub x: Vec<f64>,
    pub converged: bool,
    pub iterations: usize,
    pub reason: &'static str,
    pub residuals: Vec<f64>,
    pub estimated_residuals: Vec<f64>,
    pub iterates: Vec<Vec<f64>>,
    pub restarts: Vec<usize>,
}
impl CSRMatrix {
    pub fn gmres(&self, b: &[f64], o: GMRESOptions) -> Result<GMRESResult, MatrixError> {
        let n = self.rows();
        if self.cols() != n || b.len() != n || b.iter().any(|v| !v.is_finite()) {
            return Err(MatrixError::new(
                "GMRES requires square matrix and finite matching vector",
            ));
        }
        if o.restart < 1
            || o.restart > 1024
            || o.max_iterations > 100000
            || !o.relative_tolerance.is_finite()
            || o.relative_tolerance < 0.
            || o.relative_tolerance >= 1.
            || !o.absolute_tolerance.is_finite()
            || o.absolute_tolerance < 0.
        {
            return Err(MatrixError::new("invalid GMRES options"));
        }
        let mut diag = vec![1.; n];
        if o.jacobi {
            let rp = self.row_offsets();
            let ci = self.column_indices();
            let v = self.values();
            for i in 0..n {
                let mut found = false;
                for p in rp[i]..rp[i + 1] {
                    if ci[p] == i {
                        diag[i] = v[p];
                        found = true;
                        break;
                    }
                }
                if !found || diag[i] == 0. {
                    return Err(MatrixError::new("Jacobi requires nonzero diagonal"));
                }
            }
        }
        let norm = |v: &[f64]| v.iter().fold(0_f64, |s, &x| s.hypot(x));
        let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
        let mut x = vec![0.; n];
        let mut r = b.to_vec();
        let mut history = vec![norm(&r)];
        let mut estimates = history.clone();
        let mut frames = if o.capture { vec![x.clone()] } else { vec![] };
        let mut restarts = vec![];
        let threshold = o.absolute_tolerance.max(o.relative_tolerance * history[0]);
        let m = o.restart.min(n).min(o.max_iterations);
        let reason = 'solve: {
            if !history[0].is_finite() {
                break 'solve "nonfinite";
            }
            if history[0] <= threshold {
                break 'solve "converged";
            }
            while history.len() - 1 < o.max_iterations {
                if history.len() > 1 {
                    restarts.push(history.len() - 1);
                }
                let base = x.clone();
                let beta = norm(&r);
                let mut basis = vec![r.iter().map(|v| v / beta).collect::<Vec<f64>>()];
                let mut h = vec![vec![0.; m]; m + 1];
                let mut cs = vec![0.; m];
                let mut sn = vec![0.; m];
                let mut g = vec![0.; m + 1];
                g[0] = beta;
                let steps = m.min(o.max_iterations - (history.len() - 1));
                for j in 0..steps {
                    let z: Vec<f64> = basis[j]
                        .iter()
                        .enumerate()
                        .map(|(i, v)| v / diag[i])
                        .collect();
                    let mut w = match self.matvec(&z) {
                        Ok(v) => v,
                        Err(_) => break 'solve "nonfinite",
                    };
                    let original = norm(&w);
                    for _ in 0..2 {
                        for k in 0..=j {
                            let dot = basis[k].iter().zip(&w).map(|(a, b)| a * b).sum::<f64>();
                            h[k][j] += dot;
                            for i in 0..n {
                                w[i] -= dot * basis[k][i];
                            }
                        }
                    }
                    let tail = norm(&w);
                    if !original.is_finite()
                        || !tail.is_finite()
                        || (0..=j).any(|k| !h[k][j].is_finite())
                    {
                        break 'solve "nonfinite";
                    }
                    let happy = tail <= 8. * f64::EPSILON * original;
                    h[j + 1][j] = if happy { 0. } else { tail };
                    if !happy {
                        basis.push(w.iter().map(|v| v / tail).collect());
                    }
                    for k in 0..j {
                        let top = cs[k] * h[k][j] + sn[k] * h[k + 1][j];
                        h[k + 1][j] = -sn[k] * h[k][j] + cs[k] * h[k + 1][j];
                        h[k][j] = top;
                    }
                    let pivot = h[j][j].hypot(h[j + 1][j]);
                    if !pivot.is_finite() {
                        break 'solve "nonfinite";
                    }
                    if pivot == 0. {
                        break 'solve "breakdown";
                    }
                    cs[j] = h[j][j] / pivot;
                    sn[j] = h[j + 1][j] / pivot;
                    h[j][j] = pivot;
                    h[j + 1][j] = 0.;
                    g[j + 1] = -sn[j] * g[j];
                    g[j] = cs[j] * g[j];
                    let mut y = g[..j + 1].to_vec();
                    for k in (0..=j).rev() {
                        if h[k][k] == 0. {
                            break 'solve "breakdown";
                        }
                        let sum = (k + 1..=j).map(|q| h[k][q] * y[q]).sum::<f64>();
                        y[k] = (y[k] - sum) / h[k][k];
                    }
                    let candidate: Vec<f64> = (0..n)
                        .map(|i| {
                            base[i] + (0..=j).map(|k| basis[k][i] * y[k]).sum::<f64>() / diag[i]
                        })
                        .collect();
                    if !finite(&y) || !finite(&candidate) || !g[j + 1].is_finite() {
                        break 'solve "nonfinite";
                    }
                    let ax = match self.matvec(&candidate) {
                        Ok(v) => v,
                        Err(_) => break 'solve "nonfinite",
                    };
                    let res: Vec<f64> = b.iter().zip(ax).map(|(v, a)| v - a).collect();
                    let length = norm(&res);
                    if !length.is_finite() {
                        break 'solve "nonfinite";
                    }
                    x = candidate;
                    r = res;
                    history.push(length);
                    estimates.push(g[j + 1].abs());
                    if o.capture {
                        frames.push(x.clone());
                    }
                    if length <= threshold {
                        break 'solve "converged";
                    }
                    if happy {
                        break 'solve "breakdown";
                    }
                }
                if x == base {
                    break 'solve "stagnation";
                }
            }
            "iteration_limit"
        };
        Ok(GMRESResult {
            x,
            converged: reason == "converged",
            iterations: history.len() - 1,
            reason,
            residuals: history,
            estimated_residuals: estimates,
            iterates: frames,
            restarts,
        })
    }
}
