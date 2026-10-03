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
    println!("solution: 1, 2"); Ok(())
}
