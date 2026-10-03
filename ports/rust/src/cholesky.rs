use crate::{CSRMatrix, MatrixError};
use std::collections::BTreeMap;
/// Owned symbolic pattern. Fill steps use -1 for original entries.
#[derive(Clone, Debug)]
pub struct SparseCholeskySymbolic {
    n: usize,
    source_rp: Vec<usize>,
    source_ci: Vec<usize>,
    rp: Vec<usize>,
    ci: Vec<usize>,
    steps: Vec<isize>,
}
#[derive(Clone, Debug)]
pub struct SparseCholesky {
    lower: CSRMatrix,
}
impl SparseCholeskySymbolic {
    pub fn new(a: &CSRMatrix) -> Result<Self, MatrixError> {
        if a.rows() != a.cols() {
            return Err(MatrixError::new("Cholesky requires square matrix"));
        }
        let n = a.rows();
        let mut g = vec![BTreeMap::new(); n];
        for i in 0..n {
            for p in a.row_offsets()[i]..a.row_offsets()[i + 1] {
                let j = a.column_indices()[p];
                if i != j {
                    g[i].insert(j, -1);
                    g[j].insert(i, -1);
                }
            }
        }
        for k in 0..n {
            let ns: Vec<usize> = g[k].keys().copied().filter(|&j| j > k).collect();
            for (u, &i) in ns.iter().enumerate() {
                for &j in &ns[..u] {
                    if !g[i].contains_key(&j) {
                        g[i].insert(j, k as isize);
                        g[j].insert(i, k as isize);
                    }
                }
            }
        }
        let mut s = Self {
            n,
            source_rp: a.row_offsets().to_vec(),
            source_ci: a.column_indices().to_vec(),
            rp: vec![0],
            ci: vec![],
            steps: vec![],
        };
        for (i, row) in g.iter().enumerate() {
            for (&j, &step) in row {
                if j < i {
                    s.ci.push(j);
                    s.steps.push(step);
                }
            }
            s.ci.push(i);
            s.steps.push(-1);
            s.rp.push(s.ci.len());
        }
        Ok(s)
    }
    pub fn size(&self) -> usize {
        self.n
    }
    pub fn nnz(&self) -> usize {
        self.ci.len()
    }
    pub fn fill_count(&self) -> usize {
        self.steps.iter().filter(|&&k| k >= 0).count()
    }
    pub fn row_offsets(&self) -> &[usize] {
        &self.rp
    }
    pub fn column_indices(&self) -> &[usize] {
        &self.ci
    }
    pub fn fill_steps(&self) -> &[isize] {
        &self.steps
    }
    pub fn factorize(&self, a: &CSRMatrix) -> Result<SparseCholesky, MatrixError> {
        if a.rows() != self.n
            || a.cols() != self.n
            || a.row_offsets() != self.source_rp
            || a.column_indices() != self.source_ci
        {
            return Err(MatrixError::new("Cholesky symbolic pattern mismatch"));
        }
        let mut rows = vec![BTreeMap::new(); self.n];
        for (i, row) in rows.iter_mut().enumerate() {
            for p in a.row_offsets()[i]..a.row_offsets()[i + 1] {
                row.insert(a.column_indices()[p], a.values()[p]);
            }
        }
        for (i, row) in rows.iter().enumerate() {
            for (&j, &v) in row {
                if v != *rows[j].get(&i).unwrap_or(&0.) {
                    return Err(MatrixError::new("Cholesky requires symmetric values"));
                }
            }
        }
        let (rp, ci) = (&self.rp, &self.ci);
        let mut v = vec![0.; ci.len()];
        for i in 0..self.n {
            for p in rp[i]..rp[i + 1] {
                let j = ci[p];
                let mut t = *rows[i].get(&j).unwrap_or(&0.);
                let (mut u, mut w) = (rp[i], rp[j]);
                while u < p && w < rp[j + 1] - 1 {
                    if ci[u] == ci[w] {
                        t -= v[u] * v[w];
                        u += 1;
                        w += 1;
                    } else if ci[u] < ci[w] {
                        u += 1;
                    } else {
                        w += 1;
                    }
                }
                if !t.is_finite() {
                    return Err(MatrixError::new("nonfinite Cholesky factor"));
                }
                if i == j {
                    if t <= 0. {
                        return Err(MatrixError::new("nonpositive Cholesky pivot"));
                    }
                    v[p] = t.sqrt();
                } else {
                    v[p] = t / v[rp[j + 1] - 1];
                }
                if !v[p].is_finite() {
                    return Err(MatrixError::new("nonfinite Cholesky factor"));
                }
            }
        }
        Ok(SparseCholesky {
            lower: CSRMatrix::new(self.n, self.n, rp, ci, &v)?,
        })
    }
}
impl SparseCholesky {
    pub fn size(&self) -> usize {
        self.lower.rows()
    }
    pub fn nnz(&self) -> usize {
        self.lower.nnz()
    }
    pub fn lower(&self) -> CSRMatrix {
        self.lower.clone()
    }
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>, MatrixError> {
        if b.len() != self.size() || b.iter().any(|z| !z.is_finite()) {
            return Err(MatrixError::new("invalid Cholesky right-hand side"));
        }
        let mut x = b.to_vec();
        let (rp, ci, v) = (
            self.lower.row_offsets(),
            self.lower.column_indices(),
            self.lower.values(),
        );
        for i in 0..self.size() {
            for p in rp[i]..rp[i + 1] - 1 {
                x[i] -= v[p] * x[ci[p]];
            }
            x[i] /= v[rp[i + 1] - 1];
            if !x[i].is_finite() {
                return Err(MatrixError::new("nonfinite Cholesky solve"));
            }
        }
        for i in (0..self.size()).rev() {
            x[i] /= v[rp[i + 1] - 1];
            if !x[i].is_finite() {
                return Err(MatrixError::new("nonfinite Cholesky solve"));
            }
            for p in rp[i]..rp[i + 1] - 1 {
                x[ci[p]] -= v[p] * x[i];
                if !x[ci[p]].is_finite() {
                    return Err(MatrixError::new("nonfinite Cholesky solve"));
                }
            }
        }
        Ok(x)
    }
}
