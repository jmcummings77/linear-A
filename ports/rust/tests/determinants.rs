use linear_a::{DeterminantAlgorithm, Matrix};

fn permutation_determinant(values: &[f64], size: usize) -> i64 {
    fn visit(values: &[f64], size: usize, permutation: &mut Vec<usize>, result: &mut i64) {
        if permutation.len() == size {
            let mut product = 1_i64;
            let mut inversions = 0;
            for i in 0..size {
                product *= values[i * size + permutation[i]] as i64;
                for j in i + 1..size {
                    if permutation[i] > permutation[j] {
                        inversions += 1;
                    }
                }
            }
            *result += if inversions % 2 == 0 {
                product
            } else {
                -product
            };
            return;
        }
        for col in 0..size {
            if !permutation.contains(&col) {
                permutation.push(col);
                visit(values, size, permutation, result);
                permutation.pop();
            }
        }
    }
    let mut result = 0;
    visit(values, size, &mut Vec::new(), &mut result);
    result
}

#[test]
fn algorithms_match_independent_seeded_permutation_oracle() {
    let mut state = 781_u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((state >> 32) % 7) as f64 - 3.0
    };
    for size in 0..=6 {
        for _ in 0..6 {
            let mut values: Vec<_> = (0..size * size).map(|_| next()).collect();
            let matrix = Matrix::new(size, size, &values).unwrap();
            let expected = permutation_determinant(&values, size) as f64;
            for algorithm in [DeterminantAlgorithm::Auto, DeterminantAlgorithm::Lu] {
                let got = matrix.determinant_with(algorithm).unwrap();
                assert!(
                    (got - expected).abs() <= 1e-10 * expected.abs().max(1.0),
                    "size{size} {algorithm:?}: {got} vs{expected}"
                );
            }
            assert_eq!(matrix.values(), values);
            for row in 0..size {
                for col in row..size {
                    let value = next() + if row == col { (4 * size) as f64 } else { 0.0 };
                    values[row * size + col] = value;
                    values[col * size + row] = value;
                }
            }
            let matrix = Matrix::new(size, size, &values).unwrap();
            let expected = permutation_determinant(&values, size) as f64;
            let got = matrix
                .determinant_with(DeterminantAlgorithm::Cholesky)
                .unwrap();
            assert!((got - expected).abs() <= 1e-10 * expected.abs().max(1.0));
            assert_eq!(matrix.values(), values);
        }
    }
}

#[test]
fn extreme_scales_and_underflowing_factors_preserve_representable_results() {
    let cases: Vec<(usize, Vec<f64>, f64)> = vec![
        (
            4,
            vec![
                1e300, 0.0, 0.0, 0.0, 0.0, 1e300, 0.0, 0.0, 0.0, 0.0, 1e-300, 0.0, 0.0, 0.0, 0.0,
                1e-300,
            ],
            1.0,
        ),
        (2, vec![1e308, 1e308, 1e-308, 2e-308], 1.0),
        (2, vec![1e-200, 1e-200, 3e-124, 6e-124], f64::from_bits(1)),
        (
            3,
            vec![
                1e308, 1e308, 1e-308, 1e-308, 2e-308, 1e308, 0.0, 0.0, 1e-308,
            ],
            1e-308,
        ),
        (
            3,
            vec![f64::MAX, 0.0, 0.0, 0.0, f64::MAX, 0.0, 0.0, 0.0, 0.0],
            0.0,
        ),
    ];
    for (size, values, expected) in cases {
        let matrix = Matrix::new(size, size, &values).unwrap();
        for algorithm in [DeterminantAlgorithm::Auto, DeterminantAlgorithm::Lu] {
            let got = matrix.determinant_with(algorithm).unwrap();
            if expected == 0.0 {
                assert_eq!(got, 0.0);
            } else {
                assert!(
                    (got / expected - 1.0).abs() <= 1e-12,
                    "{values:?} {algorithm:?}: {got} vs{expected}"
                );
            }
        }
    }
    for values in [
        vec![4.0, 2.0, 2.0, 3.0],
        vec![1e300, 0.0, 0.0, 1e-300],
        vec![f64::from_bits(1)],
    ] {
        let size = if values.len() == 1 { 1 } else { 2 };
        let matrix = Matrix::new(size, size, &values).unwrap();
        let expected = matrix.determinant().unwrap();
        let got = matrix
            .determinant_with(DeterminantAlgorithm::Cholesky)
            .unwrap();
        assert!((got / expected - 1.0).abs() <= 1e-12);
    }
}

#[test]
fn cholesky_rejects_nonsymmetric_singular_and_indefinite_input() {
    for values in [
        [4.0, 2.0, 1.0, 3.0],
        [1.0, 2.0, 2.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [-1.0, 0.0, 0.0, -1.0],
    ] {
        assert!(Matrix::new(2, 2, &values)
            .unwrap()
            .determinant_with(DeterminantAlgorithm::Cholesky)
            .is_err());
    }
    assert!(Matrix::zeros(2, 3)
        .unwrap()
        .determinant_with(DeterminantAlgorithm::Cholesky)
        .is_err());
}
