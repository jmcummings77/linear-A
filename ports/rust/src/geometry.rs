use crate::{finite, Matrix, MatrixError};

fn vector3(a: &Matrix) -> bool {
    (a.rows == 3 && a.cols == 1) || (a.rows == 1 && a.cols == 3)
}

impl Matrix {
    /// Right-handed cross product; operands are 3x1 or 1x3 vectors. Retains the
    /// left shape and leaves both inputs unchanged. Nonfinite arithmetic fails.
    pub fn cross(&self, other: &Matrix) -> Result<Self, MatrixError> {
        if !vector3(self) || !vector3(other) {
            return Err(MatrixError::new(
                "cross product requires 3x1 or 1x3 vectors",
            ));
        }
        let a = &self.values;
        let b = &other.values;
        let mut result = Self::zeros(self.rows, self.cols)?;
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            result.values[i] = finite(a[j] * b[k] - a[k] * b[j])?;
        }
        Ok(result)
    }

    /// Counterclockwise active rotation of 2D column vectors, in radians.
    pub fn rotation_2d(radians: f64) -> Result<Self, MatrixError> {
        finite(radians)?;
        let (s, c) = radians.sin_cos();
        Self::new(2, 2, &[c, -s, s, c])
    }
    /// Right-handed active rotation about the positive x axis, in radians.
    pub fn rotation_x(radians: f64) -> Result<Self, MatrixError> {
        rodrigues(1.0, 0.0, 0.0, radians)
    }
    /// Right-handed active rotation about the positive y axis, in radians.
    pub fn rotation_y(radians: f64) -> Result<Self, MatrixError> {
        rodrigues(0.0, 1.0, 0.0, radians)
    }
    /// Right-handed active rotation about the positive z axis, in radians.
    pub fn rotation_z(radians: f64) -> Result<Self, MatrixError> {
        rodrigues(0.0, 0.0, 1.0, radians)
    }
    /// Right-handed active column-vector rotation about a finite nonzero 3D axis.
    /// Accepts either vector orientation and normalizes without changing it.
    pub fn rotation_axis_angle(axis: &Matrix, radians: f64) -> Result<Self, MatrixError> {
        if !vector3(axis) {
            return Err(MatrixError::new(
                "rotation axis requires a 3x1 or 1x3 vector",
            ));
        }
        let scale = axis.values.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        if scale == 0.0 {
            return Err(MatrixError::new("rotation axis must be nonzero"));
        }
        let (x, y, z) = (
            axis.values[0] / scale,
            axis.values[1] / scale,
            axis.values[2] / scale,
        );
        let norm = (x * x + y * y + z * z).sqrt();
        rodrigues(x / norm, y / norm, z / norm, radians)
    }
}

fn rodrigues(x: f64, y: f64, z: f64, radians: f64) -> Result<Matrix, MatrixError> {
    finite(radians)?;
    let (s, c) = radians.sin_cos();
    let half = (radians / 2.0).sin();
    let t = if radians.abs() < 1.0 {
        2.0 * half * half
    } else {
        1.0 - c
    };
    Matrix::new(
        3,
        3,
        &[
            c + x * x * t,
            x * y * t - z * s,
            x * z * t + y * s,
            y * x * t + z * s,
            c + y * y * t,
            y * z * t - x * s,
            z * x * t - y * s,
            z * y * t + x * s,
            c + z * z * t,
        ],
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(a: &Matrix, b: &Matrix) {
        assert_eq!((a.rows, a.cols), (b.rows, b.cols));
        for (x, y) in a.values.iter().zip(&b.values) {
            assert!((x - y).abs() < 3e-14, "{x} != {y}");
        }
    }
    fn axis(v: &[f64]) -> Matrix {
        Matrix::new(3, 1, v).unwrap()
    }
    #[test]
    fn cross_shapes_handedness_and_errors() {
        for row_a in [false, true] {
            for row_b in [false, true] {
                let a = Matrix::new(
                    if row_a { 1 } else { 3 },
                    if row_a { 3 } else { 1 },
                    &[1.0, 2.0, 3.0],
                )
                .unwrap();
                let b = Matrix::new(
                    if row_b { 1 } else { 3 },
                    if row_b { 3 } else { 1 },
                    &[4.0, 5.0, 6.0],
                )
                .unwrap();
                let result = a.cross(&b).unwrap();
                assert_eq!((result.rows, result.cols), (a.rows, a.cols));
                assert_eq!(result.values, vec![-3.0, 6.0, -3.0]);
                assert_eq!(a.values, vec![1.0, 2.0, 3.0]);
            }
        }
        assert_eq!(
            axis(&[1.0, 0.0, 0.0])
                .cross(&axis(&[0.0, 1.0, 0.0]))
                .unwrap()
                .values,
            vec![0.0, 0.0, 1.0]
        );
        assert!(axis(&[1.0, 0.0, 0.0])
            .cross(&Matrix::identity(3).unwrap())
            .is_err());
        assert!(axis(&[f64::MAX, 0.0, 0.0])
            .cross(&axis(&[0.0, 2.0, 0.0]))
            .is_err());
    }
    #[test]
    fn rotations_are_right_handed_and_orthogonal() {
        let quarter = std::f64::consts::FRAC_PI_2;
        near(
            &Matrix::rotation_x(quarter)
                .unwrap()
                .multiply(&axis(&[0.0, 1.0, 0.0]))
                .unwrap(),
            &axis(&[0.0, 0.0, 1.0]),
        );
        near(
            &Matrix::rotation_y(quarter)
                .unwrap()
                .multiply(&axis(&[0.0, 0.0, 1.0]))
                .unwrap(),
            &axis(&[1.0, 0.0, 0.0]),
        );
        near(
            &Matrix::rotation_z(quarter)
                .unwrap()
                .multiply(&axis(&[1.0, 0.0, 0.0]))
                .unwrap(),
            &axis(&[0.0, 1.0, 0.0]),
        );
        for angle in [0.0, 1e-12, -0.5, 3.14, f64::MAX] {
            for create in [
                Matrix::rotation_2d,
                Matrix::rotation_x,
                Matrix::rotation_y,
                Matrix::rotation_z,
            ] {
                let r = create(angle).unwrap();
                near(
                    &r.transpose().unwrap().multiply(&r).unwrap(),
                    &Matrix::identity(r.rows).unwrap(),
                );
                near(&create(-angle).unwrap(), &r.transpose().unwrap());
                assert!((r.determinant().unwrap() - 1.0).abs() < 3e-14);
            }
        }
    }
    #[test]
    fn axis_scaling_invariance_and_invalid_inputs() {
        let v = axis(&[1.0, 2.0, 3.0]);
        let expected = Matrix::rotation_axis_angle(&v, 0.5).unwrap();
        for scale in [1.0, 1e300, 1e-300, f64::from_bits(1)] {
            let a = axis(&[scale, 2.0 * scale, 3.0 * scale]);
            for input in [a.clone(), a.transpose().unwrap()] {
                let r = Matrix::rotation_axis_angle(&input, 0.5).unwrap();
                near(&r, &expected);
                near(&r.multiply(&v).unwrap(), &v);
                near(
                    &Matrix::rotation_axis_angle(&input, -0.5)
                        .unwrap()
                        .multiply(&r)
                        .unwrap(),
                    &Matrix::identity(3).unwrap(),
                );
                near(
                    &Matrix::rotation_axis_angle(&input, 0.0).unwrap(),
                    &Matrix::identity(3).unwrap(),
                );
            }
        }
        assert!(Matrix::rotation_axis_angle(&axis(&[0.0; 3]), 0.0).is_err());
        assert!(Matrix::rotation_axis_angle(&Matrix::identity(3).unwrap(), 0.5).is_err());
        for angle in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(Matrix::rotation_2d(angle).is_err());
            assert!(Matrix::rotation_axis_angle(&v, angle).is_err());
        }
    }
}
