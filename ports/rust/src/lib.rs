//! Owned, row-major `f64` matrices with checked shapes and arithmetic.

use std::error::Error;
use std::fmt;

mod determinant;
pub use determinant::DeterminantAlgorithm;
mod eigen_general;
pub use eigen_general::GeneralEigenDecomposition;
mod eigen;
pub use eigen::SymmetricEigenDecomposition;
mod geometry;
mod solve;
pub use solve::{FactorAlgorithm, Factorization};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MatrixError(String);

impl MatrixError {
    fn new(message: &str) -> Self {
        Self(message.to_owned())
    }
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl Error for MatrixError {}

fn finite(value: f64) -> Result<f64, MatrixError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MatrixError::new(
            "matrix arithmetic produced a nonfinite value",
        ))
    }
}

fn count(rows: usize, cols: usize) -> Result<usize, MatrixError> {
    rows.checked_mul(cols)
        .filter(|&n| n <= isize::MAX as usize / std::mem::size_of::<f64>())
        .ok_or_else(|| MatrixError::new("matrix dimensions overflow storage size"))
}

fn allocate(length: usize) -> Result<Vec<f64>, MatrixError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|_| MatrixError::new("unable to allocate matrix storage"))?;
    values.resize(length, 0.0);
    Ok(values)
}

/// Matrix values are private; construction copies its input and all operations
/// return independent storage. Empty rectangular shapes are supported.
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    values: Vec<f64>,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize, values: &[f64]) -> Result<Self, MatrixError> {
        let length = count(rows, cols)?;
        if values.len() != length {
            return Err(MatrixError::new(
                "value count does not match matrix dimensions",
            ));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(MatrixError::new("matrix values must be finite"));
        }
        let mut storage = allocate(length)?;
        storage.copy_from_slice(values);
        Ok(Self {
            rows,
            cols,
            values: storage,
        })
    }

    pub fn zeros(rows: usize, cols: usize) -> Result<Self, MatrixError> {
        Ok(Self {
            rows,
            cols,
            values: allocate(count(rows, cols)?)?,
        })
    }

    pub fn identity(size: usize) -> Result<Self, MatrixError> {
        let mut result = Self::zeros(size, size)?;
        for index in 0..size {
            result.values[index * size + index] = 1.0;
        }
        Ok(result)
    }

    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    fn offset(&self, row: usize, col: usize) -> Result<usize, MatrixError> {
        if row >= self.rows || col >= self.cols {
            return Err(MatrixError::new("matrix index is out of range"));
        }
        Ok(row * self.cols + col)
    }

    pub fn get(&self, row: usize, col: usize) -> Result<f64, MatrixError> {
        Ok(self.values[self.offset(row, col)?])
    }

    pub fn set(&mut self, row: usize, col: usize, value: f64) -> Result<(), MatrixError> {
        let offset = self.offset(row, col)?;
        self.values[offset] = finite(value)?;
        Ok(())
    }

    /// Extract a row into independent storage, including an empty row of an Nx0 matrix.
    pub fn row(&self, row: usize) -> Result<Vec<f64>, MatrixError> {
        if row >= self.rows {
            return Err(MatrixError::new("row index is out of range"));
        }
        Ok(self.values[row * self.cols..(row + 1) * self.cols].to_vec())
    }

    pub fn column(&self, col: usize) -> Result<Vec<f64>, MatrixError> {
        if col >= self.cols {
            return Err(MatrixError::new("column index is out of range"));
        }
        Ok((0..self.rows)
            .map(|row| self.values[row * self.cols + col])
            .collect())
    }

    fn same_shape(&self, other: &Self) -> Result<(), MatrixError> {
        if self.rows != other.rows || self.cols != other.cols {
            return Err(MatrixError::new("matrix dimensions must match"));
        }
        Ok(())
    }

    pub fn add(&self, other: &Self) -> Result<Self, MatrixError> {
        self.same_shape(other)?;
        let mut result = Self::zeros(self.rows, self.cols)?;
        for (index, value) in result.values.iter_mut().enumerate() {
            *value = finite(self.values[index] + other.values[index])?;
        }
        Ok(result)
    }

    pub fn subtract(&self, other: &Self) -> Result<Self, MatrixError> {
        self.same_shape(other)?;
        let mut result = Self::zeros(self.rows, self.cols)?;
        for (index, value) in result.values.iter_mut().enumerate() {
            *value = finite(self.values[index] - other.values[index])?;
        }
        Ok(result)
    }

    pub fn scale(&self, scalar: f64) -> Result<Self, MatrixError> {
        finite(scalar)?;
        let mut result = Self::zeros(self.rows, self.cols)?;
        for (index, value) in result.values.iter_mut().enumerate() {
            *value = finite(self.values[index] * scalar)?;
        }
        Ok(result)
    }

    pub fn transpose(&self) -> Result<Self, MatrixError> {
        let mut result = Self::zeros(self.cols, self.rows)?;
        for row in 0..self.rows {
            for col in 0..self.cols {
                result.values[col * self.rows + row] = self.values[row * self.cols + col];
            }
        }
        Ok(result)
    }

    pub fn multiply(&self, other: &Self) -> Result<Self, MatrixError> {
        if self.cols != other.rows {
            return Err(MatrixError::new("left columns must equal right rows"));
        }
        let mut result = Self::zeros(self.rows, other.cols)?;
        for row in 0..self.rows {
            for inner in 0..self.cols {
                let left = self.values[row * self.cols + inner];
                for col in 0..other.cols {
                    let offset = row * other.cols + col;
                    result.values[offset] = finite(
                        result.values[offset] + left * other.values[inner * other.cols + col],
                    )?;
                }
            }
        }
        Ok(result)
    }

    fn require_square(&self) -> Result<(), MatrixError> {
        if self.rows != self.cols {
            return Err(MatrixError::new("operation requires a square matrix"));
        }
        Ok(())
    }

    pub fn trace(&self) -> Result<f64, MatrixError> {
        self.require_square()?;
        let mut result = 0.0;
        for index in 0..self.rows {
            result = finite(result + self.values[index * self.cols + index])?;
        }
        Ok(result)
    }

    /// Both classifications are true for a diagonal matrix (including 0x0).
    /// A rectangular matrix has neither classification. Zeros are exact.
    pub fn triangular(&self) -> (bool, bool) {
        if self.rows != self.cols {
            return (false, false);
        }
        let mut upper = true;
        let mut lower = true;
        for row in 0..self.rows {
            for col in 0..self.cols {
                if row > col && self.values[row * self.cols + col] != 0.0 {
                    upper = false;
                }
                if row < col && self.values[row * self.cols + col] != 0.0 {
                    lower = false;
                }
            }
        }
        (upper, lower)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: usize, cols: usize, values: &[f64]) -> Matrix {
        Matrix::new(rows, cols, values).unwrap()
    }

    #[test]
    fn arithmetic_and_rectangular_product() {
        let a = matrix(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let b = matrix(3, 2, &[7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);
        assert_eq!(
            a.multiply(&b).unwrap().values(),
            &[58.0, 64.0, 139.0, 154.0]
        );
        assert_eq!(a.add(&a).unwrap(), a.scale(2.0).unwrap());
        assert_eq!(a.subtract(&a).unwrap(), Matrix::zeros(2, 3).unwrap());
        assert_eq!(
            a.transpose().unwrap().values(),
            &[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]
        );
        assert_eq!(a.transpose().unwrap().transpose().unwrap(), a);
        assert_eq!(a.multiply(&Matrix::identity(3).unwrap()).unwrap(), a);
    }

    #[test]
    fn copying_and_indexing_do_not_alias() {
        let mut source = vec![1.0, 2.0, 3.0, 4.0];
        let a = Matrix::new(2, 2, &source).unwrap();
        source[0] = 99.0;
        let mut copied = a.clone();
        copied.set(0, 0, 42.0).unwrap();
        assert_eq!(a.get(0, 0).unwrap(), 1.0);
        assert_eq!(a.row(1).unwrap(), vec![3.0, 4.0]);
        assert_eq!(a.column(1).unwrap(), vec![2.0, 4.0]);
        assert_eq!(a.values(), &[1.0, 2.0, 3.0, 4.0]);
        assert!(a.get(2, 0).is_err());
        assert!(copied.set(0, 2, 1.0).is_err());
        assert!(a.row(2).is_err());
        assert!(a.column(2).is_err());
    }

    #[test]
    fn determinant_pivoting_singularity_and_sign() {
        let swapped = matrix(2, 2, &[0.0, 2.0, 3.0, 4.0]);
        assert_eq!(swapped.determinant().unwrap(), -6.0);
        assert_eq!(swapped.trace().unwrap(), 4.0);
        let a = matrix(3, 3, &[6.0, 1.0, 1.0, 4.0, -2.0, 5.0, 2.0, 8.0, 7.0]);
        assert!((a.determinant().unwrap() + 306.0).abs() < 1e-12);
        assert_eq!(
            matrix(2, 2, &[1.0, 2.0, 2.0, 4.0]).determinant().unwrap(),
            0.0
        );
        assert_eq!(matrix(1, 1, &[-7.0]).determinant().unwrap(), -7.0);
        assert_eq!(
            matrix(2, 2, &[1e-200, 0.0, 0.0, 1.0])
                .determinant()
                .unwrap(),
            1e-200
        );
        assert_eq!(a.values(), &[6.0, 1.0, 1.0, 4.0, -2.0, 5.0, 2.0, 8.0, 7.0]);
    }

    #[test]
    fn empty_shapes_and_triangular_classification() {
        let empty = Matrix::zeros(0, 0).unwrap();
        assert_eq!(empty.trace().unwrap(), 0.0);
        assert_eq!(empty.determinant().unwrap(), 1.0);
        assert_eq!(empty.triangular(), (true, true));
        assert_eq!(Matrix::identity(3).unwrap().triangular(), (true, true));
        assert_eq!(
            matrix(2, 2, &[1.0, 2.0, 0.0, 3.0]).triangular(),
            (true, false)
        );
        assert_eq!(
            matrix(2, 2, &[1.0, 0.0, 2.0, 3.0]).triangular(),
            (false, true)
        );
        assert_eq!(Matrix::zeros(2, 0).unwrap().triangular(), (false, false));
        let product = Matrix::zeros(2, 0)
            .unwrap()
            .multiply(&Matrix::zeros(0, 3).unwrap())
            .unwrap();
        assert_eq!(product, Matrix::zeros(2, 3).unwrap());
        assert_eq!(
            Matrix::zeros(2, 0).unwrap().row(1).unwrap(),
            Vec::<f64>::new()
        );
        assert_eq!(
            Matrix::zeros(0, 2).unwrap().column(1).unwrap(),
            Vec::<f64>::new()
        );
    }

    #[test]
    fn invalid_shapes_values_and_overflow_are_errors() {
        assert!(Matrix::new(2, 2, &[1.0]).is_err());
        assert!(Matrix::zeros(usize::MAX, 2).is_err());
        assert!(Matrix::new(1, 1, &[f64::NAN]).is_err());
        let rectangular = Matrix::zeros(2, 3).unwrap();
        let square = Matrix::identity(2).unwrap();
        assert!(rectangular.add(&square).is_err());
        assert!(rectangular.subtract(&square).is_err());
        assert!(rectangular.multiply(&square).is_err());
        assert!(rectangular.trace().is_err());
        assert!(rectangular.determinant().is_err());
        assert!(square.scale(f64::INFINITY).is_err());
        let huge = matrix(1, 1, &[f64::MAX]);
        assert!(huge.add(&huge).is_err());
        assert!(huge.scale(2.0).is_err());
        assert!(huge.multiply(&huge).is_err());
        assert!(matrix(2, 2, &[f64::MAX, 0.0, 0.0, f64::MAX])
            .trace()
            .is_err());
        assert!(matrix(2, 2, &[f64::MAX, 0.0, 0.0, 2.0])
            .determinant()
            .is_err());
    }
}
