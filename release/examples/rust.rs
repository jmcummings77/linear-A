use linear_a::Matrix;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a=Matrix::new(2,2,&[4.,1.,2.,3.])?;
    let b=Matrix::new(2,1,&[6.,8.])?;
    let x=a.solve(&b)?;
    let reconstructed=a.multiply(&x)?;
    for i in 0..2 {
        assert!((x.values()[i]-(i+1) as f64).abs()<1e-12);
        assert!((reconstructed.values()[i]-b.values()[i]).abs()<1e-12);
    }
    let decomposition = a.svd()?;
    assert_eq!(decomposition.values.len(), 2);
    assert!(decomposition.values[1] > 0.);
    let inverse=a.pseudoinverse()?;let minimum=a.solve_minimum_norm(&b)?;
    assert!((inverse.get(0,0)?-0.3).abs()<1e-12 && (minimum.get(1,0)?-2.).abs()<1e-12);
    assert_eq!(a.spectral_diagnostics()?.rank,2);
    let ridge=a.solve_ridge(&b,1.)?;
    assert!((ridge.get(0,0)?-140./131.).abs()<1e-12 && (ridge.get(1,0)?-230./131.).abs()<1e-12);
    let sparse=linear_a::CSRMatrix::new(2,2,&[0,2,4],&[0,1,0,1],&[4.,1.,1.,3.])?;
    assert_eq!(sparse.matvec(&[1.,2.])?,vec![6.,7.]);
    let cg=sparse.conjugate_gradient(&[6.,7.],linear_a::CGOptions{jacobi:true,capture:true,..Default::default()})?;
    assert!(cg.converged&&(cg.x[0]-1.).abs()<1e-12&&(cg.x[1]-2.).abs()<1e-12);
    println!("solution: 1, 2"); Ok(())
}
