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
