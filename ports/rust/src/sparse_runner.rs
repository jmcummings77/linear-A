use linear_a::{CGOptions, CSRMatrix, GMRESOptions, Matrix, ILU0};
use std::io::{self, Read};
use std::time::Instant;
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
pub fn run(args: &[String]) -> Result<()> {
    if args.len()
        != (if args.get(1).map(|v| v.as_str()) == Some("gmres") {
            12
        } else {
            11
        })
    {
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
    if ![
        "spmv",
        "dense",
        "cg",
        "gmres",
        "ilu_setup",
        "ilu_apply",
        "rcm",
        "permute",
        "permutation_check",
        "rcm_solve",
        "ilu_solve",
    ]
    .contains(&op)
        || jacobi > 3
        || (op != "gmres" && jacobi > 1)
        || capture > 1
    {
        return Err("invalid sparse options".into());
    }
    let mut raw = String::new();
    io::stdin().read_to_string(&mut raw)?;
    let tokens: Vec<&str> = raw.split_whitespace().collect();
    let count = if op == "cg" || op == "gmres" {
        rows
    } else {
        cols
    };
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
    let factor = if op == "ilu_apply" || (op == "gmres" && jacobi == 2) {
        Some(ILU0::new(&a)?)
    } else {
        None
    };
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
        if op == "rcm_solve" || op == "ilu_solve" {
            let p = if op == "rcm_solve" {
                a.reverse_cuthill_mckee()?
            } else {
                (0..rows).collect()
            };
            let reordered = if op == "rcm_solve" {
                Some(a.permute_symmetric(&p)?)
            } else {
                None
            };
            let q = reordered.as_ref().unwrap_or(&a);
            let rhs = if op == "rcm_solve" {
                CSRMatrix::permute_vector(&p, &b, false)?
            } else {
                b.clone()
            };
            let options = GMRESOptions {
                restart: 20,
                relative_tolerance: rtol,
                absolute_tolerance: atol,
                max_iterations: limit,
                ..Default::default()
            };
            let result = q.gmres_preconditioned(&rhs, options, &ILU0::new(q)?)?;
            if !result.converged {
                return Err("ordering solve failed".into());
            }
            let mut out = vec![result.iterations as f64];
            out.extend(if op == "rcm_solve" {
                CSRMatrix::permute_vector(&p, &result.x, true)?
            } else {
                result.x
            });
            return Ok(out);
        }
        if op == "rcm" {
            return Ok(a
                .reverse_cuthill_mckee()?
                .iter()
                .map(|&i| i as f64)
                .collect());
        }
        if op == "permute" || op == "permutation_check" {
            if b.iter()
                .any(|x| !x.is_finite() || *x < 0. || *x >= rows as f64 || x.fract() != 0.)
            {
                return Err("invalid permutation".into());
            }
            let p: Vec<usize> = b.iter().map(|x| *x as usize).collect();
            let q = a.permute_symmetric(&p)?;
            if op == "permute" {
                let mut out: Vec<f64> = q
                    .row_offsets()
                    .iter()
                    .chain(q.column_indices().iter())
                    .map(|&x| x as f64)
                    .collect();
                out.extend(q.values());
                return Ok(out);
            }
            let x: Vec<f64> = (1..=rows).map(|i| i as f64).collect();
            let y = CSRMatrix::permute_vector(&p, &x, false)?;
            let mut out = y.clone();
            out.extend(CSRMatrix::permute_vector(&p, &y, true)?);
            out.extend(CSRMatrix::permute_vector(&p, &q.matvec(&y)?, true)?);
            return Ok(out);
        }
        if op == "ilu_setup" {
            let f = ILU0::new(&a)?;
            return Ok(vec![f.size() as f64, f.nnz() as f64]);
        }
        if op == "ilu_apply" {
            return Ok(factor.as_ref().unwrap().apply(&b)?);
        }
        if op == "spmv" {
            return Ok(a.matvec(&b)?);
        }
        if let Some((ref d, ref right)) = dense {
            return Ok(d.multiply(right)?.values().to_vec());
        }
        if op == "gmres" {
            let current = if jacobi == 3 {
                Some(ILU0::new(&a)?)
            } else {
                None
            };
            let options = GMRESOptions {
                restart: args[11].parse()?,
                relative_tolerance: rtol,
                absolute_tolerance: atol,
                max_iterations: limit,
                jacobi: jacobi == 1,
                capture: capture != 0,
            };
            let r = if let Some(p) = current.as_ref().or(factor.as_ref()) {
                a.gmres_preconditioned(&b, options, p)?
            } else {
                a.gmres(&b, options)?
            };
            if iterations > 0 {
                if !r.converged {
                    return Err("benchmark GMRES did not converge".into());
                }
                return Ok(r.x);
            }
            let reason = match r.reason {
                "converged" => 0,
                "iteration_limit" => 1,
                "breakdown" => 2,
                "nonfinite" => 3,
                _ => 4,
            };
            let mut out = vec![reason as f64, r.iterations as f64, r.residuals.len() as f64];
            out.extend(r.x);
            out.extend(r.residuals);
            out.extend(r.estimated_residuals);
            out.push(r.restarts.len() as f64);
            out.extend(r.restarts.iter().map(|&v| v as f64));
            for frame in r.iterates {
                out.extend(frame);
            }
            return Ok(out);
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
