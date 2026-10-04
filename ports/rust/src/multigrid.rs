use crate::{CSRMatrix, MatrixError};
/// A fixed SPD preconditioner. Apply preserves the RHS and returns matching length.
pub trait SymmetricPreconditioner {
    fn size(&self) -> usize;
    fn apply(&self, b: &[f64]) -> Result<Vec<f64>, MatrixError>;
}
impl SymmetricPreconditioner for crate::IC0 {
    fn size(&self) -> usize {
        self.size()
    }
    fn apply(&self, b: &[f64]) -> Result<Vec<f64>, MatrixError> {
        self.apply(b)
    }
}
/// Symmetric V-cycle for the unit five-point Dirichlet Laplacian.
#[derive(Clone, Debug)]
pub struct GeometricMultigrid {
    width: usize,
}
impl GeometricMultigrid {
    pub fn new(width: usize) -> Result<Self, MatrixError> {
        if width < 1 || width > 255 || width & (width + 1) != 0 {
            return Err(MatrixError::new("grid width must be 2^k-1 in 1..255"));
        }
        Ok(Self { width })
    }
    pub fn width(&self) -> usize {
        self.width
    }
    pub fn size(&self) -> usize {
        self.width * self.width
    }
    pub fn levels(&self) -> usize {
        (self.width + 1).trailing_zeros() as usize
    }
    pub fn matrix(&self) -> Result<CSRMatrix, MatrixError> {
        let w = self.width;
        let mut rp = vec![0];
        let mut ci = vec![];
        let mut v = vec![];
        for i in 0..self.size() {
            let mut cols = vec![];
            if i >= w {
                cols.push(i - w)
            }
            if i % w > 0 {
                cols.push(i - 1)
            }
            cols.push(i);
            if i % w + 1 < w {
                cols.push(i + 1)
            }
            if i + w < self.size() {
                cols.push(i + w)
            }
            for j in cols {
                ci.push(j);
                v.push(if j == i { 4. } else { -1. });
            }
            rp.push(v.len());
        }
        CSRMatrix::new(self.size(), self.size(), &rp, &ci, &v)
    }
    fn mv(w: usize, x: &[f64]) -> Vec<f64> {
        (0..x.len())
            .map(|i| {
                4. * x[i]
                    - if i % w > 0 { x[i - 1] } else { 0. }
                    - if i % w + 1 < w { x[i + 1] } else { 0. }
                    - if i >= w { x[i - w] } else { 0. }
                    - if i + w < x.len() { x[i + w] } else { 0. }
            })
            .collect()
    }
    fn cycle(w: usize, b: &[f64]) -> Result<Vec<f64>, MatrixError> {
        if w == 1 {
            return Ok(vec![b[0] / 4.]);
        }
        let mut x = vec![0.; b.len()];
        let smooth = |x: &mut Vec<f64>| {
            for _ in 0..2 {
                let ax = Self::mv(w, x);
                for i in 0..x.len() {
                    x[i] += (b[i] - ax[i]) / 6.;
                }
            }
        };
        smooth(&mut x);
        let ax = Self::mv(w, &x);
        let r: Vec<f64> = b.iter().zip(ax).map(|(b, a)| b - a).collect();
        let c = w / 2;
        let mut bc = vec![0.; c * c];
        for y in 0..c {
            for j in 0..c {
                for dy in 0..3 {
                    for dx in 0..3 {
                        bc[y * c + j] += (if dy == 1 { 1. } else { 0.5 })
                            * (if dx == 1 { 1. } else { 0.5 })
                            * r[(2 * y + dy) * w + 2 * j + dx];
                    }
                }
            }
        }
        let ec = Self::cycle(c, &bc)?;
        for y in 0..c {
            for j in 0..c {
                for dy in 0..3 {
                    for dx in 0..3 {
                        x[(2 * y + dy) * w + 2 * j + dx] += (if dy == 1 { 1. } else { 0.5 })
                            * (if dx == 1 { 1. } else { 0.5 })
                            * ec[y * c + j];
                    }
                }
            }
        }
        smooth(&mut x);
        if x.iter().any(|v| !v.is_finite()) {
            return Err(MatrixError::new("nonfinite multigrid cycle"));
        }
        Ok(x)
    }
    pub fn apply(&self, b: &[f64]) -> Result<Vec<f64>, MatrixError> {
        if b.len() != self.size() || b.iter().any(|v| !v.is_finite()) {
            return Err(MatrixError::new("invalid multigrid RHS"));
        }
        Self::cycle(self.width, b)
    }
}
impl SymmetricPreconditioner for GeometricMultigrid {
    fn size(&self) -> usize {
        self.size()
    }
    fn apply(&self, b: &[f64]) -> Result<Vec<f64>, MatrixError> {
        self.apply(b)
    }
}
