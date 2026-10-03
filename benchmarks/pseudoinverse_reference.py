"""Analytic inverses and minimum-norm solutions with independent Penrose checks."""
import math
from reference import matrix, assert_result


def pseudoinverse_fixtures():
    cases = []
    def add(name, op, a, expected=None, b=None, cutoff=1e-12, invalid=False, truncated=False):
        cases.append(dict(name='SVD inverse '+name, op=op, a=a, b=b, scalar=cutoff,
                          expected=expected, invalid=invalid, truncated=truncated))
    examples = [
        ('wide', matrix(1,3,[3,0,4]), matrix(3,1,[3/25,0,4/25])),
        ('tall', matrix(3,1,[3,0,4]), matrix(1,3,[3/25,0,4/25])),
        ('rank one', matrix(2,3,[1,2,3,2,4,6]), matrix(3,2,[1/70,2/70,2/70,4/70,3/70,6/70])),
        ('signed diagonal', matrix(3,3,[-4,0,0,0,0,0,0,0,2]), matrix(3,3,[-.25,0,0,0,0,0,0,0,.5])),
        ('full rank', matrix(2,2,[2,1,0,2]), matrix(2,2,[.5,-.25,0,.5])),
    ]
    for m,n in [(0,0),(0,3),(4,0),(3,2),(2,4)]:
        examples.append(('zero %dx%d'%(m,n),matrix(m,n,[0.]*(m*n)),matrix(n,m,[0.]*(m*n))))
    for name,a,p in examples:
        add(name,'pseudoinverse',a,p)
    for scale in [1e-300,1e300]:
        add('uniform scale %g'%scale,'pseudoinverse',matrix(2,2,[2*scale,0,0,scale]),matrix(2,2,[.5/scale,0,0,1/scale]))
    diagonal=matrix(2,2,[4,0,0,1])
    for cutoff, rank in [(0,2),(.24999999999999997,2),(.25,1),(.25000000000000006,1),(1,0)]:
        add('cutoff %r'%cutoff,'pseudoinverse',diagonal,matrix(2,2,[.25 if rank else 0,0,0,1 if rank==2 else 0]),cutoff=cutoff,truncated=rank<2)
        add('rank cutoff %r'%cutoff,'spectral_diagnostics',diagonal,matrix(1,3,[rank,.25,.25 if rank==2 else 1 if rank else 0]),cutoff=cutoff)
    add('retain subnormal ratio','pseudoinverse',matrix(2,2,[1,0,0,1e-300]),matrix(2,2,[1,0,0,1e300]),cutoff=0)
    add('zero spectrum','spectral_diagnostics',matrix(2,3,[0]*6),matrix(1,3,[0,0,0]))
    add('empty spectrum','spectral_diagnostics',matrix(0,3,[]),matrix(1,3,[0,0,0]))
    add('underdetermined','solve_minimum_norm',matrix(1,3,[1,1,0]),matrix(3,1,[1,1,0]),b=matrix(1,1,[2]))
    add('inconsistent rank deficiency','solve_minimum_norm',matrix(2,2,[1,1,2,2]),matrix(2,1,[.7,.7]),b=matrix(2,1,[1,3]))
    add('multiple right sides','solve_minimum_norm',matrix(1,3,[1,1,0]),matrix(3,2,[1,2,1,2,0,0]),b=matrix(1,2,[2,4]))
    add('tall least squares','solve_minimum_norm',matrix(3,2,[1,0,0,1,1,1]),matrix(2,1,[4/3,7/3]),b=matrix(3,1,[1,2,4]))
    add('empty rows','solve_minimum_norm',matrix(0,3,[]),matrix(3,2,[0]*6),b=matrix(0,2,[]))
    add('empty columns','solve_minimum_norm',matrix(3,0,[]),matrix(0,2,[]),b=matrix(3,2,[1]*6))
    add('unrepresentable inverse but finite solve','solve_minimum_norm',matrix(1,1,[1e-310]),matrix(1,1,[1]),b=matrix(1,1,[1e-310]))
    add('rhs scaling underflow','solve_minimum_norm',matrix(2,2,[1,0,0,1]),b=matrix(2,1,[1e300,1e-300]),invalid=True)
    add('inverse overflow','pseudoinverse',matrix(1,1,[1e-310]),invalid=True,cutoff=0)
    add('incompatible rhs','solve_minimum_norm',matrix(2,3,[1]*6),b=matrix(3,1,[1]*3),invalid=True)
    for cutoff in [-.1,1.1]:
        add('invalid cutoff %g'%cutoff,'pseudoinverse',diagonal,cutoff=cutoff,invalid=True)
        add('invalid diagnostics cutoff %g'%cutoff,'spectral_diagnostics',diagonal,cutoff=cutoff,invalid=True)
    return cases


def assert_inverse_result(actual, case):
    expected=case['expected']
    assert_result(actual,expected)
    # Relative checks prevent tiny expected nonzeros from passing as zero.
    for x,y in zip(actual['values'],expected['values']):
        if y: assert abs(x/y-1)<2e-9,'incorrect analytic inverse/solution'
    if case['op']=='pseudoinverse' and not case.get('truncated'):
        a=case['a'];m,n=a['rows'],a['cols']
        scale=max(map(abs,a['values']),default=0) or 1
        av=[x/scale for x in a['values']];pv=[x*scale for x in actual['values']]
        def mul(x,r,k,y,c):return [math.fsum(x[i*k+p]*y[p*c+j] for p in range(k)) for i in range(r) for j in range(c)]
        def close(x,y):
            assert max((abs(u-v) for u,v in zip(x,y)),default=0)<2e-9*max(1,max(map(abs,y),default=0)),'Penrose identity failed'
        ap=mul(av,m,n,pv,m);pa=mul(pv,n,m,av,n)
        close(mul(ap,m,m,av,n),av);close(mul(pa,n,n,pv,m),pv)
        close(ap,[ap[j*m+i] for i in range(m) for j in range(m)])
        close(pa,[pa[j*n+i] for i in range(n) for j in range(n)])
    if case['op']=='solve_minimum_norm':
        a,b=case['a'],case['b'];m,n,c=a['rows'],a['cols'],b['cols']
        # Normal-equation stationarity checked independently of the factorization.
        residual=[math.fsum(a['values'][i*n+p]*actual['values'][p*c+j] for p in range(n))-b['values'][i*c+j] for i in range(m) for j in range(c)]
        for p in range(n):
            for j in range(c):
                assert abs(math.fsum(a['values'][i*n+p]*residual[i*c+j] for i in range(m)))<2e-9,'least-squares residual not orthogonal'
