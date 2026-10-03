"""Small exact masked-Doolittle oracle, independent of the sparse elimination kernels."""
from fractions import Fraction
from sparse_reference import csr,dense_csr,multiply
from gmres_reference import transport

def reference_apply(a,b,exact=True):
    scalar=(lambda value:Fraction(str(value))) if exact else float
    n=a['rows'];pattern={(i,a['indices'][p]) for i in range(n) for p in range(a['offsets'][i],a['offsets'][i+1])}
    entries={(i,a['indices'][p]):scalar(a['values'][p]) for i in range(n) for p in range(a['offsets'][i],a['offsets'][i+1])}
    l=[[scalar(int(i==j)) for j in range(n)] for i in range(n)];u=[[scalar(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if (i,j) not in pattern:continue
            value=entries[i,j]-sum(l[i][k]*u[k][j] for k in range(min(i,j)))
            if j<i:l[i][j]=value/u[j][j]
            else:u[i][j]=value
    y=[]
    for i in range(n):y.append(scalar(b[i])-sum(l[i][j]*y[j] for j in range(i)))
    x=[scalar(0)]*n
    for i in range(n-1,-1,-1):x[i]=(y[i]-sum(u[i][j]*x[j] for j in range(i+1,n)))/u[i][i]
    return [float(v) for v in x]

def fixtures():
    cases=[]
    def add(name,a,b,op='ilu_apply',invalid=False,expected=None,reason='converged',mode=2,**options):
        if not invalid and expected is None and op!='gmres':expected=[a['rows'],len(a['values'])] if op=='ilu_setup' else reference_apply(a,b)
        cases.append(dict(name='ILU0 '+name,a=a,b=b,op=op,expected=expected,invalid=invalid,reason=reason,options={**dict(restart=3,rtol=1e-10,atol=0,limit=300,jacobi=mode if op=='gmres' else 0,capture=1 if op=='gmres' else 0),**options}))
    a=dense_csr([[4,0,1],[1,4,0],[0,1,4]])
    add('dropped fill',a,[1,2,3]);add('second RHS',a,[-3,5,2]);add('setup shape',a,[0,0,0],op='ilu_setup')
    explicit=csr(3,3,[0,2,5,7],[0,2,0,1,2,1,2],[4,1,1,4,0,1,4])
    add('explicit zero retains fill',explicit,[1,2,3])
    add('dense exact solve',dense_csr([[4,2,1],[1,3,1],[2,1,5]]),[7,5,8])
    add('negative pivots',dense_csr([[-2,1],[0,3]]),[0,6]);add('empty',csr(0,0,[0],[],[]),[])
    add('missing diagonal',dense_csr([[0,1],[1,2]]),[1,2],invalid=True)
    add('stored zero pivot',csr(1,1,[0,1],[0],[0]),[1],invalid=True)
    add('cancellation pivot',dense_csr([[1,1],[1,1]]),[2,2],invalid=True)
    add('nonsingular zero pivot',csr(2,2,[0,2,4],[0,1,0,1],[0,1,1,2]),[1,3],invalid=True)
    add('rectangular',dense_csr([[1,2]]),[1,2],invalid=True)
    add('factor overflow',dense_csr([[1e-308,1e308],[1e308,1]]),[1,1],op='ilu_setup',invalid=True)
    add('apply overflow',dense_csr([[1e-308]]),[1e308],invalid=True)
    for mode in (2,3):
        add('GMRES exact mode '+str(mode),dense_csr([[3,1],[0,2]]),[5,4],op='gmres',expected=[1,2],mode=mode,restart=1)
        for contrast in (0,3):
            a=transport(3,contrast,.2,4,30)
            add('GMRES transport %d mode %d'%(contrast,mode),a,multiply(a,[1.]*9),op='gmres',expected=[1.]*9,mode=mode)
    add('GMRES nonfinite preconditioner solve',dense_csr([[1e-320]]),[1],op='gmres',expected=None,reason='nonfinite',mode=2)
    cases[-1]['expected']=None
    return cases
