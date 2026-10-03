use linear_a::{
    DeterminantAlgorithm, GeneralEigenDecomposition, Matrix, SymmetricEigenDecomposition,
};
use std::error::Error;
use std::hint::black_box;
use std::io::{self, Read};
use std::time::Instant;

type RunnerResult<T> = Result<T, Box<dyn Error>>;

#[derive(Clone, Copy)]
enum Operation {
    Add,
    Subtract,
    Scale,
    Transpose,
    Multiply,
    Svd,
    SvdOneSweep,
    Solve,
    SolveCholesky,
    LeastSquares,
    Rcond,
    Cross,
    Rotation2D,
    Rotation3D,
    Trace,
    Determinant,
    DeterminantLu,
    DeterminantSpdLu,
    DeterminantCholesky,
    EigenSymmetric,
    EigenGeneral,
    Triangular,
}

impl Operation {
    fn parse(value: &str) -> RunnerResult<Self> {
        match value {
            "add" => Ok(Self::Add),
            "subtract" => Ok(Self::Subtract),
            "scale" => Ok(Self::Scale),
            "transpose" => Ok(Self::Transpose),
            "multiply" => Ok(Self::Multiply),
            "svd" => Ok(Self::Svd),
            "svd_one_sweep" => Ok(Self::SvdOneSweep),
            "solve" => Ok(Self::Solve),
            "solve_cholesky" => Ok(Self::SolveCholesky),
            "least_squares" => Ok(Self::LeastSquares),
            "rcond" => Ok(Self::Rcond),
            "cross" => Ok(Self::Cross),
            "rotation2d" => Ok(Self::Rotation2D),
            "rotation3d" => Ok(Self::Rotation3D),
            "trace" => Ok(Self::Trace),
            "determinant" => Ok(Self::Determinant),
            "determinant_lu" => Ok(Self::DeterminantLu),
            "determinant_spd_lu" => Ok(Self::DeterminantSpdLu),
            "determinant_cholesky" => Ok(Self::DeterminantCholesky),
            "eigen_symmetric" => Ok(Self::EigenSymmetric),
            "eigen_general" => Ok(Self::EigenGeneral),
            "triangular" => Ok(Self::Triangular),
            _ => Err(format!("unknown operation: {value}").into()),
        }
    }

    fn binary(self) -> bool {
        matches!(
            self,
            Self::Add
                | Self::Subtract
                | Self::Multiply
                | Self::Cross
                | Self::Solve
                | Self::SolveCholesky
                | Self::LeastSquares
        )
    }
}

enum Outcome {
    Matrix(Matrix),
    Scalar(f64),
    Eigen(SymmetricEigenDecomposition),
    General(GeneralEigenDecomposition),
    Triangular(bool, bool),
}

impl Outcome {
    fn checksum(&self) -> RunnerResult<f64> {
        let value = match self {
            Self::Matrix(matrix) => {
                let values = matrix.values();
                if values.is_empty() {
                    0.0
                } else {
                    values[0] + values[values.len() / 2] + values[values.len() - 1]
                }
            }
            Self::Scalar(value) => *value,
            Self::General(result) => {
                result
                    .eigenvalues_real
                    .iter()
                    .zip(&result.eigenvalues_imag)
                    .enumerate()
                    .map(|(i, (r, im))| (i + 1) as f64 * (r + im.abs()))
                    .sum::<f64>()
                    + result
                        .eigenvectors_real
                        .values()
                        .iter()
                        .zip(result.eigenvectors_imag.values())
                        .map(|(r, im)| r * r + im * im)
                        .sum::<f64>()
            }
            Self::Eigen(result) => {
                result
                    .eigenvalues
                    .iter()
                    .enumerate()
                    .map(|(i, value)| (i + 1) as f64 * value)
                    .sum::<f64>()
                    + result
                        .eigenvectors
                        .values()
                        .iter()
                        .map(|value| value * value)
                        .sum::<f64>()
            }
            Self::Triangular(..) => return Err("triangular is not a benchmark operation".into()),
        };
        if !value.is_finite() {
            return Err("checksum overflowed".into());
        }
        Ok(value)
    }

    fn print(&self) {
        match self {
            Self::Matrix(matrix) => println!(
                "{{\"rows\":{},\"cols\":{},\"values\":{:?}}}",
                matrix.rows(),
                matrix.cols(),
                matrix.values(),
            ),
            Self::Scalar(value) => println!("{{\"value\":{value}}}"),
            Self::General(result) => println!(
                "{{\"eigenvalues_real\":{:?},\"eigenvalues_imag\":{:?},\"eigenvectors_real\":{{\"rows\":{},\"cols\":{},\"values\":{:?}}},\"eigenvectors_imag\":{{\"rows\":{},\"cols\":{},\"values\":{:?}}}}}",
                result.eigenvalues_real,result.eigenvalues_imag,result.eigenvectors_real.rows(),result.eigenvectors_real.cols(),result.eigenvectors_real.values(),result.eigenvectors_imag.rows(),result.eigenvectors_imag.cols(),result.eigenvectors_imag.values()),
            Self::Eigen(result) => println!(
                "{{\"eigenvalues\":{:?},\"eigenvectors\":{{\"rows\":{},\"cols\":{},\"values\":{:?}}}}}",
                result.eigenvalues, result.eigenvectors.rows(), result.eigenvectors.cols(), result.eigenvectors.values()),
            Self::Triangular(upper, lower) => println!("{{\"upper\":{upper},\"lower\":{lower}}}"),
        }
    }
}

fn execute(
    operation: Operation,
    a: &Matrix,
    b: Option<&Matrix>,
    scalar: f64,
) -> RunnerResult<Outcome> {
    let second = || b.ok_or("binary operation requires a second matrix");
    Ok(match operation {
        Operation::Add => Outcome::Matrix(a.add(second()?)?),
        Operation::Subtract => Outcome::Matrix(a.subtract(second()?)?),
        Operation::Scale => Outcome::Matrix(a.scale(scalar)?),
        Operation::Transpose => Outcome::Matrix(a.transpose()?),
        Operation::Multiply => Outcome::Matrix(a.multiply(second()?)?),
        Operation::Svd | Operation::SvdOneSweep => {
            let r = a.svd_with(
                1e-12,
                if matches!(operation, Operation::Svd) {
                    100
                } else {
                    1
                },
            )?;
            let mut values = r.u.values().to_vec();
            values.extend(&r.values);
            values.extend(r.vt.transpose()?.values());
            Outcome::Matrix(Matrix::new(
                a.rows() + 1 + a.cols(),
                a.rows().min(a.cols()),
                &values,
            )?)
        }
        Operation::Solve => Outcome::Matrix(a.solve(second()?)?),
        Operation::SolveCholesky => Outcome::Matrix(a.factor_cholesky()?.solve(second()?)?),
        Operation::LeastSquares => Outcome::Matrix(a.least_squares(second()?)?),
        Operation::Rcond => Outcome::Scalar(a.factor_lu()?.reciprocal_condition()?),
        Operation::Cross => Outcome::Matrix(a.cross(second()?)?),
        Operation::Rotation2D => {
            if a.rows() != 0 || a.cols() != 0 {
                return Err("rotation2d requires the empty 0x0 input shape".into());
            }
            Outcome::Matrix(Matrix::rotation_2d(scalar)?)
        }
        Operation::Rotation3D => Outcome::Matrix(Matrix::rotation_axis_angle(a, scalar)?),
        Operation::Trace => Outcome::Scalar(a.trace()?),
        Operation::Determinant => Outcome::Scalar(a.determinant()?),
        Operation::DeterminantLu | Operation::DeterminantSpdLu => {
            Outcome::Scalar(a.determinant_with(DeterminantAlgorithm::Lu)?)
        }
        Operation::DeterminantCholesky => {
            Outcome::Scalar(a.determinant_with(DeterminantAlgorithm::Cholesky)?)
        }
        Operation::EigenSymmetric => Outcome::Eigen(a.eigen_symmetric()?),
        Operation::EigenGeneral => Outcome::General(a.eigen_general()?),
        Operation::Triangular => {
            let (upper, lower) = a.triangular();
            Outcome::Triangular(upper, lower)
        }
    })
}

fn integer(value: &str, name: &str) -> RunnerResult<usize> {
    if value.is_empty() || !value.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(format!("{name} must be a nonnegative integer").into());
    }
    value
        .parse()
        .map_err(|_| format!("{name} is too large").into())
}

fn number(value: &str) -> RunnerResult<f64> {
    let result: f64 = value
        .parse()
        .map_err(|_| format!("invalid number: {value}"))?;
    if !result.is_finite() {
        return Err("numbers must be finite".into());
    }
    Ok(result)
}

fn length(rows: usize, cols: usize) -> RunnerResult<usize> {
    rows.checked_mul(cols)
        .filter(|&n| n <= isize::MAX as usize / std::mem::size_of::<f64>())
        .ok_or_else(|| "matrix dimensions overflow storage size".into())
}

fn check(args: &[String]) -> RunnerResult<()> {
    if args.len() < 4 {
        return Err("usage: check OP ROWS COLS [BROWS BCOLS | SCALAR]".into());
    }
    let operation = Operation::parse(&args[1])?;
    let expected = if operation.binary() {
        6
    } else if matches!(
        operation,
        Operation::Scale | Operation::Rotation2D | Operation::Rotation3D
    ) {
        5
    } else {
        4
    };
    if args.len() != expected {
        return Err("wrong number of arguments for check operation".into());
    }
    let rows = integer(&args[2], "rows")?;
    let cols = integer(&args[3], "columns")?;
    let a_length = length(rows, cols)?;
    let (b_rows, b_cols) = if operation.binary() {
        (
            integer(&args[4], "right rows")?,
            integer(&args[5], "right columns")?,
        )
    } else {
        (0, 0)
    };
    let total = a_length
        .checked_add(length(b_rows, b_cols)?)
        .ok_or("input count overflow")?;
    let scalar = if matches!(
        operation,
        Operation::Scale | Operation::Rotation2D | Operation::Rotation3D
    ) {
        number(&args[4])?
    } else {
        1.25
    };
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let values: Vec<f64> = input
        .split_whitespace()
        .map(number)
        .collect::<RunnerResult<_>>()?;
    if values.len() != total {
        return Err(format!("expected {total} input values, received {}", values.len()).into());
    }
    let a = Matrix::new(rows, cols, &values[..a_length])?;
    let b = if operation.binary() {
        Some(Matrix::new(b_rows, b_cols, &values[a_length..])?)
    } else {
        None
    };
    execute(operation, &a, b.as_ref(), scalar)?.print();
    Ok(())
}

fn generated(size: usize, seed: usize) -> RunnerResult<Matrix> {
    let mut result = Matrix::zeros(size, size)?;
    for row in 0..size {
        for col in 0..size {
            let index = row * size + col;
            let numerator = ((index % 101 * 17 + seed % 101 * 13) % 101) as i32 - 50;
            result.set(row, col, f64::from(numerator) / 16.0)?;
        }
    }
    Ok(result)
}

fn bench(args: &[String]) -> RunnerResult<()> {
    if args.len() != 5 {
        return Err("usage: bench OP SIZE ITERATIONS SEED".into());
    }
    let operation = Operation::parse(&args[1])?;
    if matches!(operation, Operation::Triangular) {
        return Err("triangular is not a benchmark operation".into());
    }
    let size = integer(&args[2], "size")?;
    let iterations = integer(&args[3], "iterations")?;
    let seed = integer(&args[4], "seed")?;
    if size == 0 || iterations == 0 {
        return Err("size and iterations must be positive".into());
    }
    if (matches!(operation, Operation::Cross | Operation::Rotation3D) && size != 3)
        || (matches!(operation, Operation::Rotation2D) && size != 2)
    {
        return Err("cross/rotation3d require size 3; rotation2d requires size 2".into());
    }
    if seed > 2_147_483_646 {
        return Err("seed must be between 0 and 2147483646".into());
    }
    let mut a = generated(size, seed)?;
    let mut b = generated(size, seed + 1)?;
    if matches!(operation, Operation::Cross) {
        a = Matrix::new(3, 1, &a.values()[..3])?;
        let values = b.values();
        b = Matrix::new(3, 1, &[values[0], values[1], -values[2]])?;
    } else if matches!(operation, Operation::Rotation2D) {
        a = Matrix::zeros(0, 0)?;
    } else if matches!(operation, Operation::Rotation3D) {
        a = Matrix::new(3, 1, &[1.0, 2.0, 3.0])?;
    }
    let scalar = if matches!(operation, Operation::Rotation2D | Operation::Rotation3D) {
        0.5
    } else {
        1.25
    };
    if matches!(operation, Operation::EigenGeneral) {
        for i in 0..size {
            for j in 0..size {
                let k = i / 2;
                let a0 = 1.0 + (seed % 17) as f64 / 16.0 + k as f64 / 8.0;
                let b0 = 0.5 + k as f64 / 16.0;
                let value = if i == j {
                    a0
                } else if i % 2 == 0 && j == i + 1 {
                    -b0
                } else if i % 2 == 1 && j == i - 1 {
                    b0
                } else if i < j {
                    ((i * 3 + j * 5 + seed) % 11) as f64 / 32.0 - 5.0 / 32.0
                } else {
                    0.0
                };
                a.set(i, j, value)?;
            }
        }
    }
    if matches!(operation, Operation::EigenSymmetric) {
        for row in 0..size {
            for col in 0..size {
                let value = if row == col {
                    2.0 + (seed % 17) as f64 / 16.0
                } else if row.abs_diff(col) == 1 {
                    -1.0
                } else {
                    0.0
                };
                a.set(row, col, value)?;
            }
        }
    }
    if matches!(
        operation,
        Operation::DeterminantSpdLu | Operation::DeterminantCholesky
    ) {
        for row in 0..size {
            for col in row + 1..size {
                let value = (a.get(row, col)? + a.get(col, row)?) / 2.0;
                a.set(row, col, value)?;
                a.set(col, row, value)?;
            }
        }
    }
    if matches!(
        operation,
        Operation::Determinant
            | Operation::DeterminantLu
            | Operation::DeterminantSpdLu
            | Operation::DeterminantCholesky
    ) {
        for index in 0..size {
            a.set(index, index, a.get(index, index)? + size as f64 * 4.0)?;
        }
    }
    let warmups = iterations.clamp(5, 100);
    for _ in 0..warmups {
        black_box(execute(operation, black_box(&a), Some(black_box(&b)), scalar)?.checksum()?);
    }
    let started = Instant::now();
    let mut checksum = 0.0;
    for _ in 0..iterations {
        checksum += execute(operation, black_box(&a), Some(black_box(&b)), scalar)?.checksum()?;
        if !checksum.is_finite() {
            return Err("accumulated checksum overflowed".into());
        }
    }
    let elapsed_ns = started.elapsed().as_nanos();
    println!("{{\"elapsed_ns\":{elapsed_ns},\"iterations\":{iterations},\"checksum\":{checksum}}}");
    Ok(())
}

fn run() -> RunnerResult<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("check") => check(&args),
        Some("bench") => bench(&args),
        _ => Err("usage: linear-a-rust {check|bench} OP ...".into()),
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}
