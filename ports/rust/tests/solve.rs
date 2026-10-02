use linear_a::{FactorAlgorithm, Factorization, Matrix};
#[test]
fn reusable_factor_snapshots_and_multiple_rhs() {
    for algorithm in [
        FactorAlgorithm::Lu,
        FactorAlgorithm::Cholesky,
        FactorAlgorithm::Qr,
    ] {
        let factor = {
            let a = Matrix::new(2, 2, &[4.0, 1.0, 1.0, 3.0]).unwrap();
            Factorization::new(&a, algorithm).unwrap()
        };
        let rhs = Matrix::new(2, 2, &[6.0, 5.0, 7.0, 4.0]).unwrap();
        for _ in 0..3 {
            let x = factor.solve(&rhs).unwrap();
            for (actual, expected) in x.values().iter().zip([1.0, 1.0, 2.0, 1.0]) {
                assert!((actual - expected).abs() < 1e-12);
            }
        }
        assert!((factor.reciprocal_condition().unwrap() - 0.44).abs() < 1e-12);
        assert_eq!(rhs.values(), &[6.0, 5.0, 7.0, 4.0]);
        assert!(factor.solve(&Matrix::zeros(1, 1).unwrap()).is_err());
    }
}
