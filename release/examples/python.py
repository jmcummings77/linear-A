from matrix import Matrix

a = Matrix(2, 2, [4, 1, 2, 3])
b = Matrix(2, 1, [6, 8])
x = a.solve(b)
assert max(abs(u-v) for u,v in zip(x.values, [1, 2])) < 1e-12
assert max(abs(u-v) for u,v in zip(a.multiply(x).values, b.values)) < 1e-12
print("solution: 1, 2")

r = a.svd()
assert len(r.values) == 2 and r.values[0] >= r.values[1] > 0

p=a.pseudoinverse(); minimum=a.solve_minimum_norm(b); d=a.spectral_diagnostics()
assert abs(p[0,0]-.3)<1e-12 and abs(minimum[1,0]-2)<1e-12 and d.rank==2

ridge=a.solve_ridge(b,1)
assert abs(ridge[0,0]-140/131)<1e-12 and abs(ridge[1,0]-230/131)<1e-12

from matrix import CSRMatrix
sparse=CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,1,1,3])
assert sparse.matvec([1,2])==[6,7]
cg=sparse.conjugate_gradient([6,7],jacobi=True,capture=True)
assert cg.converged and max(abs(x-y) for x,y in zip(cg.x,[1,2]))<1e-12

gm=sparse.gmres([6,7],restart=2,jacobi=True,capture=True)
assert gm.converged and max(abs(x-y) for x,y in zip(gm.x,[1,2]))<1e-12

from matrix import ILU0
f=ILU0(sparse)
for rhs,expected in [([6,7],[1,2]),([11,13],[20/11,41/11])]:
    assert max(abs(x-y) for x,y in zip(f.apply(rhs),expected))<1e-12
    assert sparse.gmres(rhs,preconditioner=f).converged

p=sparse.reverse_cuthill_mckee();q=sparse.permute_symmetric(p)
y=CSRMatrix.permute_vector(p,[1,2])
assert CSRMatrix.permute_vector(p,q.matvec(y),True)==[6,7]

from matrix import SparseCholeskySymbolic
assert sparse.approximate_minimum_degree()==[0,1]
plan=SparseCholeskySymbolic(sparse)
chol=plan.factorize(sparse)
assert abs(chol.solve([6,7])[0]-1)<1e-12 and abs(chol.solve([11,13])[0]-20/11)<1e-12
assert plan.fill_count==0 and chol.lower.nnz==3

from matrix import IC0
ic=IC0(sparse)
for rhs in ([6,7],[11,13]):
    pcg=sparse.conjugate_gradient(rhs,preconditioner=ic)
    assert pcg.converged and pcg.iterations==1
