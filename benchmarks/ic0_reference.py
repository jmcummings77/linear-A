"""Dense rational LDLᵀ oracle for zero-fill incomplete Cholesky.

Uses a dense mask and exact fractions, independently of the CSR square-root kernels.
"""
from fractions import Fraction as F
import math
from sparse_reference import csr,dense_csr,multiply

BREAKDOWN=[[1.5,1,0,1],[1,1.5,1,0],[0,1,1.5,-1],[1,0,-1,1.5]]

def reference(a,exact=True):
    scalar=(lambda v:F(str(v))) if exact else float
    n=a['rows']
    if a['cols']!=n:raise ValueError('square required')
    entries=[[scalar(0) for _ in range(n)] for _ in range(n)];mask=[[i==j for j in range(n)] for i in range(n)]
    for i in range(n):
        for p in range(a['offsets'][i],a['offsets'][i+1]):
            j=a['indices'][p];entries[i][j]=scalar(a['values'][p]);mask[i][j]=mask[j][i]=True
    if any(entries[i][j]!=entries[j][i] for i in range(n) for j in range(n)):raise ValueError('symmetry required')
    l=[[scalar(int(i==j)) for j in range(n)] for i in range(n)];d=[scalar(0)]*n
    for i in range(n):
        for j in range(i):
            if mask[i][j]:l[i][j]=(entries[i][j]-sum(l[i][k]*d[k]*l[j][k] for k in range(j)))/d[j]
        d[i]=entries[i][i]-sum(l[i][k]**2*d[k] for k in range(i))
        if d[i]<=0:raise ValueError('nonpositive IC0 pivot')
    return l,d,mask

def factor(a,exact=True):
    l,d,mask=reference(a,exact);rp=[0];ci=[];v=[]
    for i in range(a['rows']):
        for j in range(i+1):
            if mask[i][j]:ci.append(j);v.append(float(l[i][j])*math.sqrt(float(d[j])))
        rp.append(len(v))
    return rp,ci,v

def apply(a,b,exact=True):
    scalar=(lambda v:F(str(v))) if exact else float
    l,d,_=reference(a,exact);n=len(d);y=[]
    for i in range(n):y.append(scalar(b[i])-sum(l[i][j]*y[j] for j in range(i)))
    x=[scalar(0)]*n
    for i in reversed(range(n)):x[i]=y[i]/d[i]-sum(l[j][i]*x[j] for j in range(i+1,n))
    return list(map(float,x))

def fixtures():
    cases=[]
    def add(name,a,b,op='ic0_factor',invalid=False,expected=None,reason='converged',**options):
        if not invalid and expected is None and op!='cg':expected=sum(factor(a),[]) if op=='ic0_factor' else apply(a,b)
        cases.append(dict(name='IC0 '+name,a=a,b=b,op=op,invalid=invalid,expected=expected,reason=reason,options={**dict(jacobi=2 if op=='cg' else 0,capture=int(op=='cg'),limit=100,rtol=1e-10,atol=0),**options}))
    a=dense_csr([[4,1,0,1],[1,4,1,0],[0,1,4,1],[1,0,1,4]])
    add('dropped fill',a,[1,2,3,4]);add('apply',a,[1,2,3,4],op='ic0_apply');add('second RHS',a,[-3,0,5,1],op='ic0_apply')
    explicit=csr(4,4,[0,3,7,10,13],[0,1,3,0,1,2,3,1,2,3,0,2,3],[4,1,1,1,4,1,0,1,4,1,1,1,4])
    add('one-sided stored zero retains fill',explicit,[1,2,3,4])
    add('empty',csr(0,0,[0],[],[]),[]);add('empty apply',csr(0,0,[0],[],[]),[],op='ic0_apply')
    for mode in (2,3):
        add('CG mode '+str(mode),a,multiply(a,[1,2,3,4]),op='cg',expected=[1,2,3,4],jacobi=mode)
        full=dense_csr([[4,1],[1,3]])
        add('CG exact one update '+str(mode),full,[6,7],op='cg',expected=[1,2],limit=1,jacobi=mode)
    add('CG iteration limit',a,[1,2,3,4],op='cg',reason='iteration_limit',limit=1)
    add('CG apply overflow',dense_csr([[1e-320]]),[1],op='cg',reason='nonfinite')
    for name,bad in [('SPD breakdown',dense_csr(BREAKDOWN)),('zero pivot',dense_csr([[1,1],[1,1]])),('negative pivot',dense_csr([[-1]])),('missing diagonal',dense_csr([[0,1],[1,2]])),('asymmetric',dense_csr([[2,1],[0,2]])),('rectangular',dense_csr([[1,2]])),('factor overflow',dense_csr([[1e-308,1e308],[1e308,1]]))]:
        add(name,bad,[1]*bad['cols'],invalid=True)
    add('apply overflow',dense_csr([[1e-320]]),[1],op='ic0_apply',invalid=True)
    return cases
