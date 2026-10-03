"""Independent dense Cholesky and Boolean elimination references for small fixtures."""
import math
from sparse_reference import csr,dense_csr,multiply

def analyze(a):
    n=a['rows'];assert n==a['cols']
    edges=[[False]*n for _ in range(n)]
    steps=[[-1]*n for _ in range(n)]
    for i in range(n):
        edges[i][i]=True
        for p in range(a['offsets'][i],a['offsets'][i+1]):
            j=a['indices'][p];edges[i][j]=edges[j][i]=True
    for k in range(n):
        for i in range(k+1,n):
            for j in range(k+1,i):
                if edges[i][k] and edges[j][k] and not edges[i][j]:
                    edges[i][j]=edges[j][i]=True;steps[i][j]=steps[j][i]=k
    rp=[0];ci=[];birth=[]
    for i in range(n):
        for j in range(i+1):
            if edges[i][j]:ci.append(j);birth.append(steps[i][j])
        rp.append(len(ci))
    return rp,ci,birth

def factor(a):
    n=a['rows'];d=[[0.]*n for _ in range(n)]
    for i in range(n):
        for p in range(a['offsets'][i],a['offsets'][i+1]):d[i][a['indices'][p]]=a['values'][p]
    l=[[0.]*n for _ in range(n)]
    for j in range(n):
        l[j][j]=math.sqrt(d[j][j]-sum(l[j][k]**2 for k in range(j)))
        for i in range(j+1,n):l[i][j]=(d[i][j]-sum(l[i][k]*l[j][k] for k in range(j)))/l[j][j]
    rp,ci,_=analyze(a)
    return rp,ci,[l[i][ci[p]] for i in range(n) for p in range(rp[i],rp[i+1])]

def fixtures():
    cases=[]
    matrices=[('empty',csr(0,0,[0],[],[])),('scalar',dense_csr([[4]])),('tiny pivot',dense_csr([[1e-280]])),('diagonal',dense_csr([[2,0,0],[0,3,0],[0,0,5]])),('path',dense_csr([[4,1,0],[1,4,1],[0,1,4]])),('star fill',dense_csr([[5,1,1,1],[1,3,0,0],[1,0,3,0],[1,0,0,3]])),('cycle',dense_csr([[4,1,0,1],[1,4,1,0],[0,1,4,1],[1,0,1,4]])),('stored zero',csr(3,3,[0,3,5,7],[0,1,2,0,1,0,2],[4,0,1,0,3,1,3])),('one-sided zero',csr(2,2,[0,2,3],[0,1,1],[4,0,3])),('numerical cancellation',dense_csr([[4,2,2],[2,5,1],[2,1,5]]))]
    for name,a in matrices:
        n=a['rows'];x=[1+(i%3)*.25 for i in range(n)];b=multiply(a,x)
        for op in ('chol_symbolic','chol_factor','chol_solve','chol_total','chol_rcm_total'):
            expected=sum(analyze(a),[]) if op=='chol_symbolic' else sum(factor(a),[]) if op=='chol_factor' else x
            cases.append(dict(name=name+' '+op,op=op,a=a,b=b,expected=expected,invalid=False))
    for name,a in [('indefinite',dense_csr([[1,2],[2,1]])),('semidefinite',dense_csr([[1,1],[1,1]])),('asymmetric',dense_csr([[4,1],[0,3]])),('missing diagonal',csr(2,2,[0,1,2],[1,0],[1,1])),('negative diagonal',dense_csr([[-1]])),('overflow',dense_csr([[1e-300,1e300],[1e300,1]]))]:
        for op in ('chol_factor','chol_solve','chol_total','chol_rcm_total'):
            cases.append(dict(name=name+' '+op,op=op,a=a,b=[1.]*a['rows'],invalid=True))
    for op in ('chol_symbolic','chol_factor','chol_solve','chol_total','chol_rcm_total'):
        cases.append(dict(name='rectangular '+op,op=op,a=csr(1,2,[0,1],[0],[1]),b=[1,1],invalid=True))
    return cases

def check_result(actual,case):
    from sparse_reference import norm
    v=actual['values'];a=case['a'];n=a['rows']
    assert set(actual)=={'rows','cols','values'} and actual['rows']==1 and actual['cols']==len(v)
    assert len(v)==len(case['expected']) and all(math.isfinite(x) for x in v)
    assert all(abs(x-y)<=1e-10*max(abs(y),1e-300) for x,y in zip(v,case['expected']))
    if case['op']=='chol_factor':
        m=(len(v)-n-1)//2;rp=v[:n+1];ci=v[n+1:n+1+m];values=v[n+1+m:]
        l=[[0.]*n for _ in range(n)];d=[[0.]*n for _ in range(n)]
        for i in range(n):
            for p in range(int(rp[i]),int(rp[i+1])):l[i][int(ci[p])]=values[p]
            for p in range(a['offsets'][i],a['offsets'][i+1]):d[i][a['indices'][p]]=a['values'][p]
        error=norm([sum(l[i][k]*l[j][k] for k in range(n))-d[i][j] for i in range(n) for j in range(n)])
        assert error<=1e-11*max(norm(a['values']),1e-300)
    elif case['op']!='chol_symbolic':
        assert norm([b-y for b,y in zip(case['b'],multiply(a,v))])<=1e-11*max(norm(case['b']),1e-300)
    return {}
