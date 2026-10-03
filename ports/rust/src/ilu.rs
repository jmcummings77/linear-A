use crate::{CSRMatrix, MatrixError};
/// Owned ILU(0). No fill, permutations or diagonal shifts.
#[derive(Clone, Debug)]
pub struct ILU0 {
    n: usize,
    rp: Vec<usize>,
    ci: Vec<usize>,
    v: Vec<f64>,
    d: Vec<usize>,
}
impl ILU0 {
    pub fn new(a: &CSRMatrix) -> Result<Self, MatrixError> {
        if a.rows() != a.cols() {
            return Err(MatrixError::new("ILU0 requires square matrix"));
        }
        let mut f = Self {
            n: a.rows(),
            rp: a.row_offsets().to_vec(),
            ci: a.column_indices().to_vec(),
            v: a.values().to_vec(),
            d: vec![0; a.rows()],
        };
        for i in 0..f.n {
            let p = f.find(i, i);
            if p == f.rp[i + 1] || f.ci[p] != i {
                return Err(MatrixError::new("ILU0 requires stored diagonal"));
            }
            f.d[i] = p;
        }
        for i in 0..f.n {
            for p in f.rp[i]..f.d[i] {
                let j = f.ci[p];
                f.v[p] /= f.v[f.d[j]];
                if !f.v[p].is_finite() {
                    return Err(MatrixError::new("nonfinite ILU0 factor"));
                }
                for q in f.d[j] + 1..f.rp[j + 1] {
                    let k = f.find(i, f.ci[q]);
                    if k < f.rp[i + 1] && f.ci[k] == f.ci[q] {
                        f.v[k] -= f.v[p] * f.v[q];
                        if !f.v[k].is_finite() {
                            return Err(MatrixError::new("nonfinite ILU0 factor"));
                        }
                    }
                }
            }
            if f.v[f.d[i]] == 0. {
                return Err(MatrixError::new("zero ILU0 pivot"));
            }
        }
        Ok(f)
    }
    fn find(&self, i: usize, j: usize) -> usize {
        let mut lo = self.rp[i];
        let mut hi = self.rp[i + 1];
        while lo < hi {
            let m = lo + (hi - lo) / 2;
            if self.ci[m] < j {
                lo = m + 1
            } else {
                hi = m
            }
        }
        lo
    }
    pub fn size(&self) -> usize {
        self.n
    }
    pub fn nnz(&self) -> usize {
        self.v.len()
    }
    pub fn apply(&self, b: &[f64]) -> Result<Vec<f64>, MatrixError> {
        if b.len() != self.n || b.iter().any(|v| !v.is_finite()) {
            return Err(MatrixError::new("invalid ILU0 vector"));
        }
        let mut x = b.to_vec();
        for i in 0..self.n {
            for p in self.rp[i]..self.d[i] {
                x[i] -= self.v[p] * x[self.ci[p]];
            }
            if !x[i].is_finite() {
                return Err(MatrixError::new("nonfinite ILU0 solve"));
            }
        }
        for i in (0..self.n).rev() {
            for p in self.d[i] + 1..self.rp[i + 1] {
                x[i] -= self.v[p] * x[self.ci[p]];
            }
            x[i] /= self.v[self.d[i]];
            if !x[i].is_finite() {
                return Err(MatrixError::new("nonfinite ILU0 solve"));
            }
        }
        Ok(x)
    }
}
