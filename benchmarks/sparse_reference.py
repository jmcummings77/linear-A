"""Shared CSR and CG fixtures, independent residual checks, and diffusion inputs."""
import math


def csr(rows,cols,rp,ci,values):return dict(rows=rows,cols=cols,offsets=rp,indices=ci,values=values)
def dense_csr(values,rows=None,cols=None):
    rows=len(values) if rows is None else rows;cols=len(values[0]) if rows and cols is None else cols or 0
    rp,ci,v=[0],[],[]
    for row in values:
        for j,x in enumerate(row):
            if x:ci.append(j);v.append(x)
        rp.append(len(v))
    return csr(rows,cols,rp,ci,v)

def multiply(a,x):return [sum(a['values'][p]*x[a['indices'][p]] for p in range(a['offsets'][i],a['offsets'][i+1])) for i in range(a['rows'])]
def norm(v):
    s=0.
    for x in v:s=math.hypot(s,x)
    return s

def fixtures():
    cases=[]
    def add(name,a,b,op='cg',expected=None,reason='converged',invalid=False,**options):
        cases.append(dict(name=name,a=a,b=b,op=op,expected=expected,reason=reason,invalid=invalid,options=dict(rtol=1e-10,atol=0,limit=100,jacobi=0,capture=1,**options) if not set(options)&{'rtol','atol','limit','jacobi','capture'} else {**dict(rtol=1e-10,atol=0,limit=100,jacobi=0,capture=1),**options}))
    a=dense_csr([[4,1],[1,3]])
    add('SPD analytic',a,[6,7],expected=[1,2]);add('SPD Jacobi',a,[6,7],expected=[1,2],jacobi=1)
    diagonal=dense_csr([[.01,0,0],[0,1,0],[0,0,100]])
    add('Jacobi diagonal one update',diagonal,[.01,2,300],expected=[1,2,3],jacobi=1,limit=1)
    add('iteration cap',a,[6,7],reason='iteration_limit',limit=0)
    add('one iteration cap',a,[6,7],reason='iteration_limit',limit=1)
    add('absolute tolerance',a,[1e-12,0],expected=[0,0],atol=1e-11)
    add('zero RHS',a,[0,0],expected=[0,0])
    add('empty CG',csr(0,0,[0],[],[]),[],expected=[])
    add('no frames',a,[6,7],expected=[1,2],capture=0)
    add('negative curvature',dense_csr([[-1,0],[0,1]]),[1,0],reason='breakdown')
    add('singular incompatible',dense_csr([[1,0],[0,0]]),[0,1],reason='breakdown')
    add('nonfinite arithmetic',dense_csr([[1]]),[1e200],reason='nonfinite')
    add('underflow breakdown',dense_csr([[1]]),[1e-200],reason='breakdown')
    for name,mat,b in [('nonsymmetric',dense_csr([[1,1],[0,1]]),[1,2]),('rectangular CG',csr(1,2,[0,1],[0],[1]),[1])]:add(name,mat,b,invalid=True)
    add('Jacobi missing diagonal',dense_csr([[0,1],[1,1]]),[1,2],jacobi=1,invalid=True)
    add('Jacobi negative diagonal',dense_csr([[-1]]),[1],jacobi=1,invalid=True)
    for key,val in [('rtol',-1),('rtol',1),('atol',-1),('limit',100001)]:add('invalid '+key+str(val),a,[6,7],invalid=True,**{key:val})
    rectangular=csr(3,4,[0,2,2,4],[0,3,1,2],[2,-1,3,4])
    add('rectangular empty row',rectangular,[1,2,3,4],op='spmv',expected=[-2,0,18])
    add('explicit zero',csr(1,2,[0,2],[0,1],[0,2]),[5,3],op='spmv',expected=[6])
    add('empty rows',csr(0,3,[0],[],[]),[1,2,3],op='spmv',expected=[])
    add('empty columns',csr(3,0,[0,0,0,0],[],[]),[],op='spmv',expected=[0,0,0])
    add('dense comparison',rectangular,[1,2,3,4],op='dense',expected=[-2,0,18])
    for name,mat in [('first offset',csr(1,1,[1,1],[0],[1])),('last offset',csr(1,1,[0,0],[0],[1])),('nonmonotone',csr(2,2,[0,2,1],[0],[1])),('duplicate',csr(1,2,[0,2],[0,0],[1,2])),('unsorted',csr(1,2,[0,2],[1,0],[1,2])),('column out of range',csr(1,1,[0,1],[1],[1])),('negative column',csr(1,1,[0,1],[-1],[1])),('negative offset',csr(1,1,[0,-1],[],[]))]:
        add('invalid CSR '+name,mat,[1]*mat['cols'],op='spmv',invalid=True)
    add('nonfinite product',csr(1,1,[0,1],[0],[1e308]),[2],op='spmv',invalid=True)
    for contrast in [0,2]:
        mat=diffusion(4,contrast);truth=[math.sin(i+1) for i in range(16)];rhs=multiply(mat,truth)
        for jacobi in [0,1]:add('diffusion contrast %s Jacobi %s'%(contrast,jacobi),mat,rhs,expected=truth,jacobi=jacobi,limit=200)
    from graph_reference import fixtures as graph_fixtures
    cases.extend(graph_fixtures())
    return cases

def protocol(case,iterations=0):
    a=case['a'];o={**dict(rtol=1e-10,atol=0,limit=1000,jacobi=0,capture=0),**case.get('options',{})}
    args=['sparse',case['op'],str(a['rows']),str(a['cols']),str(len(a['values'])),str(iterations),str(o['rtol']),str(o['atol']),str(o['limit']),str(o['jacobi']),str(o['capture'])]
    data=a['offsets']+a['indices']+a['values']+case['b']
    if case['op']=='gmres':args.append(str(o.get('restart',30)))
    return args,' '.join(str(x) for x in data)

def check_result(actual,case):
    assert set(actual)=={'rows','cols','values'} and actual['rows']==1
    v=actual['values'];assert actual['cols']==len(v) and all(isinstance(x,(float,int)) and math.isfinite(x) for x in v)
    if case['op']!='cg':
        assert len(v)==len(case['expected'])
        assert all(abs(x-y)<=1e-12*max(1,abs(y)) for x,y in zip(v,case['expected']));return
    n=case['a']['rows'];reason,iterations,length=v[:3]
    assert all(int(x)==x for x in [reason,iterations,length]);reason,iterations,length=map(int,[reason,iterations,length])
    assert ['converged','iteration_limit','breakdown','nonfinite'][reason]==case['reason']
    assert length==iterations+1 and 0<=iterations<=case['options']['limit']
    x=v[3:3+n];history=v[3+n:3+n+length];flat=v[3+n+length:]
    assert len(x)==n and len(history)==length
    assert len(flat)==(n*length if case['options']['capture'] else 0)
    if case['expected'] is not None:assert len(x)==len(case['expected']) and all(abs(a-b)<1e-7*max(1,abs(b)) for a,b in zip(x,case['expected']))
    threshold=max(case['options']['atol'],case['options']['rtol']*norm(case['b']))
    if reason==0:assert norm([b-y for b,y in zip(case['b'],multiply(case['a'],x))])<=threshold*(1+1e-6)+1e-14
    if case['options']['capture']:
        for k in range(length):
            frame=flat[k*n:(k+1)*n];res=norm([b-y for b,y in zip(case['b'],multiply(case['a'],frame))])
            assert abs(history[k]-res)<=1e-11*max(res,history[0],1e-300),'incorrect true residual history'
        assert flat[-n:]==x if n else not flat


def diffusion(size,contrast=0):
    """Five-point positive diffusion with zero Dirichlet boundary; no dense allocation."""
    n=size*size;conductivity=[10**(contrast*(i%size)/max(1,size-1)) for i in range(n)]
    rp,ci,v=[0],[],[]
    for i in range(n):
        y,x=divmod(i,size);entries={};diagonal=0.
        for dy,dx in [(-1,0),(0,-1),(0,1),(1,0)]:
            yy,xx=y+dy,x+dx
            if 0<=yy<size and 0<=xx<size:
                j=yy*size+xx;weight=math.sqrt(conductivity[i]*conductivity[j]);entries[j]=-weight
            else:weight=conductivity[i]
            diagonal+=weight
        entries[i]=diagonal
        for j,value in sorted(entries.items()):ci.append(j);v.append(value)
        rp.append(len(v))
    return csr(n,n,rp,ci,v)
