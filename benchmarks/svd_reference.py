"""Analytic spectra and invariant checks; no numerical SVD used as an oracle."""
import math
import random
from reference import matrix


def svd_fixtures():
    cases=[]
    def add(name,m,n,values,spectrum=None,invalid=False,op='svd'):
        cases.append(dict(name='SVD '+name,op=op,a=matrix(m,n,values),b=None,invalid=invalid,spectrum=spectrum))
    for m,n in [(0,0),(0,3),(4,0),(3,2),(2,4)]:
        add('zero %dx%d'%(m,n),m,n,[0.]*(m*n),[0.]*min(m,n))
    add('signed diagonal',3,3,[-3,0,0,0,0,0,0,0,5],[5,3,0])
    add('repeated values',2,3,[0,2,0,-2,0,0],[2,2])
    add('row vector',1,3,[3,0,4],[5])
    add('column vector',3,1,[3,0,4],[5])
    add('rank one',4,3,[float((i+1)*(j+1)) for i in range(4) for j in range(3)],[math.sqrt(420),0,0])
    add('flat image',24,24,[1.]*576,[24]+[0]*23)
    # Rotated left/right orthonormal factors, exact rational coefficients and known spectrum.
    for scale in [1,1e-200,1e200,1e-300,1e300]:
        add('rotated spectrum %g'%scale,3,2,[scale*x for x in [2.16,-2.88,2.88,-3.84,1.6,1.2]],[6*scale,2*scale])
    add('small isolated value',2,2,[1,0,0,1e-200],[1,1e-200])
    add('small correlated column',2,2,[1,1e-150,0,1e-150],[1,1e-150])
    rng=random.Random(28491)
    for m,n in [(3,3),(5,3),(3,5),(9,6)]:
        add('generated %dx%d'%(m,n),m,n,[rng.uniform(-1,1) for _ in range(m*n)])
    add('singular-value overflow',1,2,[1.7e308,1.7e308],invalid=True)
    add('input scaling underflow',2,2,[1e300,0,0,1e-300],invalid=True)
    add('iteration exhaustion',3,3,[1,2,3,4,5,7,6,8,9],invalid=True,op='svd_one_sweep')
    return cases


def assert_svd_result(result,a,spectrum=None):
    m,n=a['rows'],a['cols'];k=min(m,n)
    assert result['rows']==m+1+n and result['cols']==k,'incorrect packed SVD shape'
    data=result['values']
    assert len(data)==(m+1+n)*k and all(type(x) in (int,float) and math.isfinite(x) for x in data),'nonfinite or malformed factors'
    u=data[:m*k];s=data[m*k:(m+1)*k];v=data[(m+1)*k:]
    assert all(x>=0 for x in s) and all(s[i]>=s[i+1] for i in range(k-1)),'singular values are not descending/nonnegative'
    for factors,count in [(u,m),(v,n)]:
        for p in range(k):
            for q in range(k):
                dot=math.fsum(factors[i*k+p]*factors[i*k+q] for i in range(count))
                assert abs(dot-(p==q))<2e-9,'singular vectors are not orthonormal'
    scale=max(map(abs,a['values']),default=0) or 1
    residual=[]
    for i in range(m):
        for j in range(n):
            reconstructed=math.fsum(u[i*k+p]*(s[p]/scale)*v[j*k+p] for p in range(k))
            residual.append(reconstructed-a['values'][i*n+j]/scale)
    assert math.sqrt(math.fsum(x*x for x in residual))<2e-9*max(1,math.sqrt(m*n)),'SVD reconstruction residual'
    if spectrum is not None:
        largest=max(spectrum,default=0)
        for actual,expected in zip(s,spectrum):
            # Nonzero analytic values use a relative check, including tiny isolated values.
            if expected: assert abs(actual/expected-1)<2e-9,'incorrect analytic singular value'
            else: assert actual<=2e-12*max(largest,1e-300),'spurious nonzero singular value'
