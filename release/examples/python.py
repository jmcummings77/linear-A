from matrix import Matrix

a = Matrix(2, 2, [4, 1, 2, 3])
b = Matrix(2, 1, [6, 8])
x = a.solve(b)
assert max(abs(u-v) for u,v in zip(x.values, [1, 2])) < 1e-12
assert max(abs(u-v) for u,v in zip(a.multiply(x).values, b.values)) < 1e-12
print("solution: 1, 2")
