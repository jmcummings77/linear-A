use super::{finite, Matrix, MatrixError};

/// Selects a floating-point determinant method.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DeterminantAlgorithm {
    /// Exact-zero triangular shortcut, otherwise partial-pivot LU.
    Auto,
    Lu,
    /// Requires exactly symmetric positive-definite input.
    Cholesky,
}

fn parts(mut value: f64) -> (f64, i32) {
    let mut bits = value.to_bits();
    let mut encoded = ((bits >> 52) & 0x7ff) as i32;
    let mut exponent = 0;
    if encoded == 0 {
        value *= 18014398509481984.0;
        bits = value.to_bits();
        encoded = ((bits >> 52) & 0x7ff) as i32;
        exponent -= 54;
    }
    exponent += encoded - 1023;
    (
        f64::from_bits((bits & 0x800fffffffffffff) | (1023_u64 << 52)),
        exponent,
    )
}

fn scale_normalized(mantissa: f64, exponent: i64) -> f64 {
    if exponent > 1023 {
        return f64::INFINITY.copysign(mantissa);
    }
    if exponent < -1075 {
        return 0.0_f64.copysign(mantissa);
    }
    if exponent < -1022 {
        return mantissa * 2.0_f64.powi((exponent + 1022) as i32) * f64::MIN_POSITIVE;
    }
    mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
}

fn scale(value: f64, exponent: i32) -> f64 {
    if value == 0.0 {
        return value;
    }
    let (mantissa, original) = parts(value);
    scale_normalized(mantissa, i64::from(original) + i64::from(exponent))
}

struct Product {
    mantissa: f64,
    exponent: i64,
    zero: bool,
}

impl Product {
    fn new() -> Self {
        Self {
            mantissa: 1.0,
            exponent: 0,
            zero: false,
        }
    }
    fn multiply(&mut self, value: f64) {
        if value == 0.0 {
            self.zero = true;
            return;
        }
        let (fraction, exponent) = parts(value);
        self.mantissa *= fraction;
        self.exponent += i64::from(exponent);
        if self.mantissa.abs() >= 2.0 {
            self.mantissa *= 0.5;
            self.exponent += 1;
        }
    }
    fn result(self) -> Result<f64, MatrixError> {
        if self.zero {
            return Ok(0.0);
        }
        finite(scale_normalized(self.mantissa, self.exponent))
    }
}

fn quotient_product(numerator: f64, pivot: f64, value: f64) -> f64 {
    if value == 0.0 {
        return 0.0;
    }
    let (a, ae) = parts(numerator);
    let (b, be) = parts(pivot);
    let (c, ce) = parts(value);
    let (fraction, exponent) = parts(a / b * c);
    scale_normalized(fraction, i64::from(ae - be + ce + exponent))
}

impl Matrix {
    /// Auto determinant; det(0x0)=1. The input is not changed.
    pub fn determinant(&self) -> Result<f64, MatrixError> {
        self.determinant_with(DeterminantAlgorithm::Auto)
    }

    /// Cubic-time elimination with quadratic working storage, except triangular
    /// Auto inputs. Reversible binary row scaling and a separate determinant
    /// exponent reduce avoidable range loss; rounding and conditioning remain.
    pub fn determinant_with(&self, algorithm: DeterminantAlgorithm) -> Result<f64, MatrixError> {
        self.require_square()?;
        for &value in &self.values {
            finite(value)?;
        }
        if algorithm == DeterminantAlgorithm::Cholesky {
            return self.cholesky_determinant();
        }
        if algorithm == DeterminantAlgorithm::Auto {
            let (upper, lower) = self.triangular();
            if upper || lower {
                let mut product = Product::new();
                for i in 0..self.rows {
                    product.multiply(self.values[i * self.cols + i]);
                }
                return product.result();
            }
        }
        let size = self.rows;
        let mut work = self.values.clone();
        let mut product = Product::new();
        for row in 0..size {
            let mut largest = 0.0_f64;
            for col in 0..size {
                largest = largest.max(work[row * size + col].abs());
            }
            if largest == 0.0 {
                return Ok(0.0);
            }
            let (_, exponent) = parts(largest);
            let preserves = (0..size).all(|col| {
                let value = work[row * size + col];
                scale(scale(value, -exponent), exponent) == value
            });
            if preserves {
                for col in 0..size {
                    work[row * size + col] = scale(work[row * size + col], -exponent);
                }
                product.exponent += i64::from(exponent);
            }
        }
        for col in 0..size {
            let mut pivot_row = col;
            for row in col + 1..size {
                if work[row * size + col].abs() > work[pivot_row * size + col].abs() {
                    pivot_row = row;
                }
            }
            if work[pivot_row * size + col] == 0.0 {
                return Ok(0.0);
            }
            if pivot_row != col {
                for j in col..size {
                    work.swap(col * size + j, pivot_row * size + j);
                }
                product.mantissa = -product.mantissa;
            }
            let pivot = work[col * size + col];
            product.multiply(pivot);
            for row in col + 1..size {
                let numerator = work[row * size + col];
                let factor = numerator / pivot;
                let scaled_update = numerator != 0.0 && factor.abs() < f64::MIN_POSITIVE;
                work[row * size + col] = 0.0;
                for j in col + 1..size {
                    let term = if scaled_update {
                        quotient_product(numerator, pivot, work[col * size + j])
                    } else {
                        factor * work[col * size + j]
                    };
                    work[row * size + j] = finite(work[row * size + j] - term)?;
                }
            }
        }
        product.result()
    }

    fn cholesky_determinant(&self) -> Result<f64, MatrixError> {
        let size = self.rows;
        let mut work = self.values.clone();
        for row in 0..size {
            for col in 0..row {
                if work[row * size + col] != work[col * size + row] {
                    return Err(MatrixError::new("Cholesky requires exact symmetry"));
                }
            }
        }
        let mut product = Product::new();
        for row in 0..size {
            for col in 0..=row {
                let mut sum = work[row * size + col];
                for k in 0..col {
                    sum -= work[row * size + k] * work[col * size + k];
                }
                if !sum.is_finite() || row == col && sum <= 0.0 {
                    return Err(MatrixError::new("Cholesky requires finite positive pivots"));
                }
                work[row * size + col] = finite(if row == col {
                    sum.sqrt()
                } else {
                    sum / work[col * size + col]
                })?;
            }
            product.multiply(work[row * size + row]);
            product.multiply(work[row * size + row]);
        }
        product.result()
    }
}
