"""Exact rational Gauss-Jordan reference for small solver checks.

Production solvers use LU, Cholesky, or Householder QR. Least-squares references
use normal equations *in exact rational arithmetic*, never floating point.
"""
from fractions import Fraction
import random
from reference import matrix


def exact_solve(a, b):
    n, p = a['rows'], b['cols']
    if a['cols'] != n or b['rows'] != n:
        raise ValueError('incompatible exact solve shapes')
    work = [[Fraction(v) for v in a['values'][i*n:(i+1)*n] + b['values'][i*p:(i+1)*p]] for i in range(n)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if work[i][k]), None)
        if pivot is None:
            raise ValueError('singular exact system')
        work[k], work[pivot] = work[pivot], work[k]
        divisor = work[k][k]
        work[k] = [v/divisor for v in work[k]]
        for i in range(n):
            if i != k:
                coefficient = work[i][k]
                work[i] = [v-coefficient*w for v, w in zip(work[i], work[k])]
    return matrix(n, p, [float(v) for row in work for v in row[n:]])


def exact_least_squares(a, b):
    m, n, p = a['rows'], a['cols'], b['cols']
    if b['rows'] != m or m < n:
        raise ValueError('incompatible least-squares shapes')
    av, bv = list(map(Fraction, a['values'])), list(map(Fraction, b['values']))
    ata = matrix(n, n, [sum(av[k*n+i]*av[k*n+j] for k in range(m)) for i in range(n) for j in range(n)])
    atb = matrix(n, p, [sum(av[k*n+i]*bv[k*p+j] for k in range(m)) for i in range(n) for j in range(p)])
    return exact_solve(ata, atb)


def solve_fixtures():
    cases = []
    def add(name, op, a, b=None, invalid=False, expected=None):
        case = dict(name=name, op=op, a=a, b=b, invalid=invalid)
        if not invalid:
            case['expected'] = expected if expected is not None else (exact_least_squares(a,b) if op == 'least_squares' else exact_solve(a,b))
        cases.append(case)
    pivot = matrix(3,3,[0,2,1,3,1,4,2,5,6])
    rhs = matrix(3,2,[1,2,3,4,5,6])
    for op in ('solve','least_squares'):
        add('solver pivoted multiple RHS '+op,op,pivot,rhs)
        add('solver empty '+op,op,matrix(0,0,[]),matrix(0,3,[]))
        add('solver zero RHS columns '+op,op,pivot,matrix(3,0,[]))
        add('solver singular '+op,op,matrix(2,2,[1,2,2,4]),matrix(2,1,[1,2]),invalid=True)
        add('solver RHS mismatch '+op,op,pivot,matrix(2,1,[1,2]),invalid=True)
        add('solver underdetermined '+op,op,matrix(1,2,[1,2]),matrix(1,1,[3]),invalid=True)
    for scale in (1,1e-300,1e300):
        a=matrix(2,2,[4*scale,scale,scale,3*scale]); b=matrix(2,2,[6*scale,5*scale,7*scale,4*scale])
        for op in ('solve','solve_cholesky','least_squares'):
            add('solver uniform scale %g %s'%(scale,op),op,a,b,expected=matrix(2,2,[1,1,2,1]))
    for name,a in [('nonsymmetric',matrix(2,2,[1,1,0,1])),('indefinite',matrix(2,2,[1,2,2,1])),('semidefinite',matrix(2,2,[1,1,1,1]))]:
        add('solver Cholesky '+name,'solve_cholesky',a,matrix(2,1,[1,1]),invalid=True)
    add('solver empty Cholesky','solve_cholesky',matrix(0,0,[]),matrix(0,2,[]))
    add('solver tall LU rejected','solve',matrix(2,1,[1,2]),matrix(2,1,[1,2]),invalid=True)
    add('solver QR empty columns','least_squares',matrix(3,0,[]),matrix(3,2,[1,2,3,4,5,6]))
    add('solver QR column pivoting and noisy observations','least_squares',matrix(4,2,[1,10,1,20,1,30,1,40]),matrix(4,2,[2,1,4,5,5,2,8,9]))
    add('solver QR numerical rank threshold','least_squares',matrix(2,2,[1,0,0,1e-18]),matrix(2,1,[1,1e-18]),invalid=True)
    add('solver LU retains small pivot','solve',matrix(2,2,[1,0,0,1e-18]),matrix(2,1,[1,1e-18]))
    add('solver unsupported mixed scaling','solve',matrix(2,2,[1e300,0,0,1e-300]),matrix(2,1,[1,1]),invalid=True)
    add('solver overflow output','solve',matrix(1,1,[1e-300]),matrix(1,1,[1e300]),invalid=True)
    add('solver subnormal scale','solve',matrix(1,1,[5e-324]),matrix(1,1,[5e-324]))
    for name,a,expected in [('diagonal',matrix(2,2,[1,0,0,4]),.25),('SPD',matrix(2,2,[4,1,1,3]),.44),('empty',matrix(0,0,[]),1),('row pivot',matrix(2,2,[0,2,1,3]),.1)]:
        add('solver rcond '+name,'rcond',a,expected={'value':expected})
    rng=random.Random(741)
    for n in range(1,6):
        a=matrix(n,n,[rng.randrange(-5,6)+(20 if i==j else 0) for i in range(n) for j in range(n)])
        b=matrix(n,3,[rng.randrange(-10,11) for _ in range(n*3)])
        for op in ('solve','least_squares'): add('solver exact seeded %d %s'%(n,op),op,a,b)
    return cases
