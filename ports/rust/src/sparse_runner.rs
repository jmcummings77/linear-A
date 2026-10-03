use linear_a::{CGOptions, CSRMatrix, Matrix};
use std::io::{self, Read};
use std::time::Instant;
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
pub fn run(args: &[String]) -> Result<()> {
    if args.len() != 11 {
        return Err("invalid sparse protocol".into());
    }
    let op = args[1].as_str();
    let rows: usize = args[2].parse()?;
    let cols: usize = args[3].parse()?;
    let nnz: usize = args[4].parse()?;
    let iterations: usize = args[5].parse()?;
    let rtol: f64 = args[6].parse()?;
    let atol: f64 = args[7].parse()?;
    let limit: usize = args[8].parse()?;
    let jacobi: usize = args[9].parse()?;
    let capture: usize = args[10].parse()?;
    if !["spmv", "dense", "cg"].contains(&op) || jacobi > 1 || capture > 1 {
        return Err("invalid sparse options".into());
    }
    let mut raw = String::new();
    io::stdin().read_to_string(&mut raw)?;
    let tokens: Vec<&str> = raw.split_whitespace().collect();
    let count = if op == "cg" { rows } else { cols };
    if tokens.len()
        != rows
            .checked_add(1)
            .and_then(|x| nnz.checked_mul(2).and_then(|v| x.checked_add(v)))
            .and_then(|x| x.checked_add(count))
            .ok_or("input size overflow")?
    {
        return Err("incorrect sparse input count".into());
    }
    let rp = tokens[..rows + 1]
        .iter()
        .map(|s| s.parse())
        .collect::<std::result::Result<Vec<usize>, _>>()?;
    let ci = tokens[rows + 1..rows + 1 + nnz]
        .iter()
        .map(|s| s.parse())
        .collect::<std::result::Result<Vec<usize>, _>>()?;
    let v = tokens[rows + 1 + nnz..rows + 1 + 2 * nnz]
        .iter()
        .map(|s| s.parse())
        .collect::<std::result::Result<Vec<f64>, _>>()?;
    let b = tokens[rows + 1 + 2 * nnz..]
        .iter()
        .map(|s| s.parse())
        .collect::<std::result::Result<Vec<f64>, _>>()?;
    let a = CSRMatrix::new(rows, cols, &rp, &ci, &v)?;
    let dense = if op == "dense" {
        let mut data = vec![0.; rows.checked_mul(cols).ok_or("size overflow")?];
        for i in 0..rows {
            for p in rp[i]..rp[i + 1] {
                data[i * cols + ci[p]] = v[p];
            }
        }
        Some((Matrix::new(rows, cols, &data)?, Matrix::new(cols, 1, &b)?))
    } else {
        None
    };
    let compute = || -> Result<Vec<f64>> {
        if op == "spmv" {
            return Ok(a.matvec(&b)?);
        }
        if let Some((ref d, ref right)) = dense {
            return Ok(d.multiply(right)?.values().to_vec());
        }
        let r = a.conjugate_gradient(
            &b,
            CGOptions {
                relative_tolerance: rtol,
                absolute_tolerance: atol,
                max_iterations: limit,
                jacobi: jacobi != 0,
                capture: capture != 0,
            },
        )?;
        if iterations > 0 {
            if !r.converged {
                return Err("benchmark CG did not converge".into());
            }
            return Ok(r.x);
        }
        let reason = match r.reason {
            "converged" => 0,
            "iteration_limit" => 1,
            "breakdown" => 2,
            _ => 3,
        };
        let mut out = vec![reason as f64, r.iterations as f64, r.residuals.len() as f64];
        out.extend(r.x);
        out.extend(r.residuals);
        for frame in r.iterates {
            out.extend(frame);
        }
        Ok(out)
    };
    if iterations == 0 {
        let values = compute()?;
        println!(
            "{{\"rows\":1,\"cols\":{},\"values\":{:?}}}",
            values.len(),
            values
        );
        return Ok(());
    }
    for _ in 0..3 {
        std::hint::black_box(compute()?);
    }
    let mut checksum = 0.;
    let start = Instant::now();
    for _ in 0..iterations {
        for x in std::hint::black_box(compute()?) {
            checksum += x;
        }
    }
    let elapsed = start.elapsed().as_nanos();
    if !checksum.is_finite() {
        return Err("nonfinite checksum".into());
    }
    println!(
        "{{\"elapsed_ns\":{},\"iterations\":{},\"checksum\":{}}}",
        elapsed, iterations, checksum
    );
    Ok(())
}
