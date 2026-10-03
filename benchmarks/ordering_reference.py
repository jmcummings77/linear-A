"""Independent ordering fixtures and permutation identities (zero-based new-to-old)."""
from sparse_reference import csr,dense_csr,multiply

def rcm(a):
    n=a['rows'];assert n==a['cols']
    # Build edge pairs, then use explicit breadth-first frontiers.
    edges=set()
    for i in range(n):
        for k in range(a['offsets'][i],a['offsets'][i+1]):
            j=a['indices'][k]
            if i!=j:edges.add((i,j));edges.add((j,i))
    graph=[{j for u,j in edges if u==i} for i in range(n)]
    key=lambda i:(len(graph[i]),i)
    remaining=set(range(n));order=[]
    while remaining:
        front=[min(remaining,key=key)];remaining.remove(front[0])
        while front:
            order+=front;following=[]
            for i in front:
                neighbors=sorted(graph[i]&remaining,key=key)
                following+=neighbors;remaining.difference_update(neighbors)
            front=following
    return order[::-1]

def permute(a,p):
    assert len(p)==a['rows']==a['cols'] and sorted(p)==list(range(len(p)))
    inverse={old:new for new,old in enumerate(p)};rp=[0];ci=[];values=[]
    for old in p:
        row={inverse[a['indices'][k]]:a['values'][k] for k in range(a['offsets'][old],a['offsets'][old+1])}
        for j in sorted(row):ci.append(j);values.append(row[j])
        rp.append(len(ci))
    return csr(len(p),len(p),rp,ci,values)

def bandwidth(a):
    return max((abs(i-a['indices'][k]) for i in range(a['rows']) for k in range(a['offsets'][i],a['offsets'][i+1])),default=0)

def fixtures():
    cases=[]
    matrices=[('empty',csr(0,0,[0],[],[])),('isolated',csr(4,4,[0,0,0,0,0],[],[])),('diagonal',dense_csr([[1,0,0],[0,2,0],[0,0,3]])),('nonsymmetric path',dense_csr([[4,0,0,2],[0,5,3,0],[0,0,6,0],[0,7,0,8]])),('explicit zero',csr(3,3,[0,2,3,4],[0,2,1,2],[4,0,5,6])),('cycle ties',dense_csr([[3,1,0,1],[1,3,1,0],[0,1,3,1],[1,0,1,3]])),('disconnected',dense_csr([[2,1,0,0,0],[0,2,0,0,0],[0,0,2,0,1],[0,0,0,2,0],[0,0,1,0,2]]))]
    for name,a in matrices:
        p=rcm(a);q=permute(a,p);n=a['rows'];x=list(range(1,n+1))
        for op,b,expected in [('rcm',[0]*n,p),('permute',p,q['offsets']+q['indices']+q['values']),('permutation_check',p,[x[i] for i in p]+x+multiply(a,x))]:
            cases.append(dict(name=name+' '+op,op=op,a=a,b=b,expected=expected,invalid=False))
    a=matrices[3][1]
    for p in ([0,0,2,3],[-1,1,2,3],[0,1,2,4],[0,.5,2,3]):
        cases.append(dict(name='invalid permutation '+str(p),op='permute',a=a,b=p,invalid=True))
    for op in ('rcm','permute'):
        cases.append(dict(name='rectangular '+op,op=op,a=csr(1,2,[0,1],[1],[2]),b=[0,1],invalid=True))
    # Not a self-inverse ordering: catches reversing the new/old convention.
    p=[2,0,3,1];q=permute(a,p);x=[1,2,3,4]
    cases.append(dict(name='non-involutory permutation',op='permutation_check',a=a,b=p,expected=[x[i] for i in p]+x+multiply(a,x),invalid=False))
    a=dense_csr([[4,1],[2,3]])
    for op in ('ilu_solve','rcm_solve'):
        cases.append(dict(name=op+' original residual',op=op,a=a,b=[6,8],expected=[1,2],invalid=False,options=dict(rtol=1e-10,atol=0,limit=100)))
    return cases

def check_solve(actual,case):
    from sparse_reference import norm
    assert set(actual)=={'rows','cols','values'} and actual['rows']==1
    v=actual['values'];n=case['a']['rows'];assert len(v)==n+1 and actual['cols']==len(v)
    assert isinstance(v[0],(int,float)) and int(v[0])==v[0] and 0<=v[0]<=case['options']['limit']
    x=v[1:];assert all(__import__('math').isfinite(t) for t in x)
    assert norm([b-y for b,y in zip(case['b'],multiply(case['a'],x))])<=max(case['options']['atol'],case['options']['rtol']*norm(case['b']))*(1+1e-5)+1e-300
    assert all(abs(a-b)<1e-7*max(1,abs(b)) for a,b in zip(x,case['expected']))
    return dict(iterations=int(v[0]))
