"""Independent analytic ridge cases; Decimal protects the extreme-scale oracle."""
from decimal import Decimal, localcontext
from reference import matrix, assert_result


def ridge_fixtures():
    cases=[]
    def add(name,a,b,lam,expected=None,invalid=False):
        cases.append(dict(name='ridge '+name,op='solve_ridge',a=a,b=b,scalar=lam,expected=expected,invalid=invalid))
    for s,b,lam in [(2,3,1),(-2,3,1),(2,3,0),(1e300,1e300,1e300),
                    (1e-300,1e300,1e300),(1e-300,1e300,1e-300),
                    (1e-300,1e-300,1e-300),(1e-160,1e-160,5e-324),
                    (1,1,1e-300),(1,1,1e300),(1e-310,1e-310,0)]:
        with localcontext() as ctx:
            ctx.prec=100
            ds,db,dl=map(Decimal.from_float,map(float,(s,b,lam)))
            expected=float(ds*db/(ds*ds+dl))
        add('scalar %g/%g/%g'%(s,b,lam),matrix(1,1,[s]),matrix(1,1,[b]),lam,matrix(1,1,[expected]))
    add('wide multiple RHS',matrix(1,3,[1,1,0]),matrix(1,2,[2,4]),2,matrix(3,2,[.5,1,.5,1,0,0]))
    add('rank deficient inconsistent',matrix(2,2,[1,1,2,2]),matrix(2,1,[1,3]),4,matrix(2,1,[.5,.5]))
    add('tall',matrix(3,2,[1,0,0,1,1,1]),matrix(3,1,[1,2,4]),1,matrix(2,1,[1.125,1.625]))
    add('zero lambda rank deficient',matrix(1,3,[1,1,0]),matrix(1,1,[2]),0,matrix(3,1,[1,1,0]))
    add('signed diagonal',matrix(2,2,[-2,0,0,1]),matrix(2,2,[3,6,4,8]),1,matrix(2,2,[-1.2,-2.4,2,4]))
    for m,n,c in [(0,3,2),(3,0,2),(3,2,0),(3,2,2),(2,3,2)]:
        add('empty/zero %d/%d/%d'%(m,n,c),matrix(m,n,[0]*(m*n)),matrix(m,c,[1]*(m*c)),1,matrix(n,c,[0]*(n*c)))
    add('negative lambda',matrix(1,1,[1]),matrix(1,1,[1]),-1,invalid=True)
    add('wrong RHS',matrix(2,3,[1]*6),matrix(1,1,[1]),1,invalid=True)
    add('RHS scale loss',matrix(2,2,[1,0,0,1]),matrix(2,1,[1e300,1e-300]),1,invalid=True)
    return cases


def assert_ridge_result(actual,case):
    assert_result(actual,case['expected'])
    for x,y in zip(actual['values'],case['expected']['values']):
        if y: assert abs(x/y-1)<2e-9,'incorrect ridge analytic solution'
    # Check Aᵀ(Ax-b)+lambda*x using high precision, independently of SVD.
    with localcontext() as ctx:
        ctx.prec=100
        a,b=case['a'],case['b'];m,n,c=a['rows'],a['cols'],b['cols']
        av,bv,x=([Decimal.from_float(float(v)) for v in vs] for vs in (a['values'],b['values'],actual['values']))
        lam=Decimal.from_float(float(case['scalar']))
        residual=[sum((av[i*n+p]*x[p*c+j] for p in range(n)),Decimal(0))-bv[i*c+j] for i in range(m) for j in range(c)]
        for p in range(n):
            for j in range(c):
                terms=[av[i*n+p]*residual[i*c+j] for i in range(m)]+[lam*x[p*c+j]]
                scale=sum((abs(av[i*n+p]*bv[i*c+j]) for i in range(m)),Decimal(0))+abs(lam*x[p*c+j])
                assert abs(sum(terms))<=Decimal('2e-9')*max(scale,Decimal('1e-630')),'ridge stationarity failed'
