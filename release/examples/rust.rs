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
    let ilu=linear_a::ILU0::new(&sparse)?;
    assert!((ilu.apply(&[6.,7.])?[0]-1.).abs()<1e-12&&(ilu.apply(&[11.,13.])?[0]-20./11.).abs()<1e-12);
    assert!(sparse.gmres_preconditioned(&[11.,13.],Default::default(),&ilu)?.converged);
    let gm=sparse.gmres(&[6.,7.],linear_a::GMRESOptions{restart:2,jacobi:true,capture:true,..Default::default()})?;
    assert!(gm.converged&&(gm.x[0]-1.).abs()<1e-12&&(gm.x[1]-2.).abs()<1e-12);
    let cg=sparse.conjugate_gradient(&[6.,7.],linear_a::CGOptions{jacobi:true,capture:true,..Default::default()})?;
    assert!(cg.converged&&(cg.x[0]-1.).abs()<1e-12&&(cg.x[1]-2.).abs()<1e-12);
    let order=sparse.reverse_cuthill_mckee()?;let reordered=sparse.permute_symmetric(&order)?;let y=linear_a::CSRMatrix::permute_vector(&order,&[1.,2.],false)?;
    assert_eq!(linear_a::CSRMatrix::permute_vector(&order,&reordered.matvec(&y)?,true)?,vec![6.,7.]);
    assert_eq!(sparse.approximate_minimum_degree()?,vec![0,1]);
    let plan=linear_a::SparseCholeskySymbolic::new(&sparse)?;let chol=plan.factorize(&sparse)?;
    assert!((chol.solve(&[6.,7.])?[0]-1.).abs()<1e-12&&(chol.solve(&[11.,13.])?[0]-20./11.).abs()<1e-12);
    assert_eq!(chol.lower().nnz(),3);
    let ic=linear_a::IC0::new(&sparse)?;
    for rhs in [[6.,7.],[11.,13.]] {let pcg=sparse.conjugate_gradient_preconditioned(&rhs,Default::default(),Some(&ic))?;assert!(pcg.converged && pcg.iterations==1);}
    let mg=linear_a::GeometricMultigrid::new(3)?;let grid=mg.matrix()?;let rhs_mg=grid.matvec(&[1.;9])?;
    assert!(grid.conjugate_gradient_preconditioned(&rhs_mg,Default::default(),Some(&mg))?.converged);
    assert_eq!(mg.apply(&rhs_mg)?.len(),9);assert_eq!(mg.levels(),2);
    println!("solution: 1, 2"); Ok(())
}
