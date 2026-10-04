"""Explicit elimination oracle and quotient-bound checks, independent of port storage."""
import random
from sparse_reference import csr,dense_csr,multiply,diffusion
from ordering_reference import permute

def graph(a):
    assert a['rows']==a['cols']
    edges=set()
    for i in range(a['rows']):
        for k in range(a['offsets'][i],a['offsets'][i+1]):
            j=a['indices'][k]
            if i!=j:edges.add(tuple(sorted((i,j))))
    return edges

def amd(a,trace=False):
    n=a['rows'];original=graph(a);edges=original.copy();active=set(range(n))
    bounds=[sum(i in e for e in edges) for i in range(n)];cliques=[];order=[];frames=[]
    while active:
        exact={i:sum(i in e and set(e)<=active for e in edges) for i in active}
        estimates={i:min(bounds[i],len(active)-1) for i in active}
        assert all(exact[i]<=estimates[i]<=len(active)-1 for i in active)
        pivot=min(active,key=lambda i:(estimates[i],i))
        neighbors={j for e in edges if pivot in e for j in e if j!=pivot and j in active}
        additions={tuple(sorted((i,j))) for i in neighbors for j in neighbors if i<j}-edges
        frames.append(dict(pivot=pivot,neighbors=sorted(neighbors),fill=sorted(additions),degree=exact[pivot],bound=estimates[pivot],estimates=estimates,exact=exact))
        edges.update(additions);active.remove(pivot);order.append(pivot)
        cliques=[c for c in cliques if pivot not in c]+[neighbors]
        # Reconstruct uncovered original edges, rather than mutate adjacency lists.
        uncovered={e for e in original if set(e)<=active and not any(set(e)<=c for c in cliques)}
        for i in neighbors:
            bounds[i]=len(neighbors)-1+sum(i in e for e in uncovered)+sum(len(c-neighbors) for c in cliques[:-1] if i in c)
    return (order,frames) if trace else order

def matrix(size,problem):
    if problem in ('grid','scrambled'):
        a=diffusion(size,2)
        if problem=='scrambled':
            p=list(range(a['rows']));random.Random(2026).shuffle(p);a=permute(a,p)
        return a
    n=size*size;edges=set()
    if problem=='tree':edges={(i,(i-1)//2) for i in range(1,n)}
    elif problem=='irregular':
        rng=random.Random(2026)
        edges={(i,rng.randrange(i)) for i in range(1,n)}
        edges.update(tuple(sorted(rng.sample(range(n),2))) for _ in range(2*n))
    else:raise ValueError('unknown problem')
    d=[[0.]*n for _ in range(n)]
    for i,j in edges:d[i][j]=d[j][i]=-1.
    for i in range(n):d[i][i]=1+sum(abs(v) for v in d[i])
    return dense_csr(d)

def fixtures():
    from ordering_reference import fixtures as ordering_cases
    cases=[]
    for c in ordering_cases():
        if c['op']=='rcm':
            c=dict(c,op='amd',name=c['name'].replace('rcm','amd'))
            if not c['invalid']:c['expected']=amd(c['a'])
            cases.append(c)
    for kind in ('grid','scrambled','tree','irregular'):
        a=matrix(4,kind);x=[1+(i%7)/10 for i in range(a['rows'])]
        cases.append(dict(name=kind+' AMD',op='amd',a=a,b=x,expected=amd(a),invalid=False))
        cases.append(dict(name=kind+' AMD solve',op='chol_amd_total',a=a,b=multiply(a,x),expected=x,invalid=False))
    # A strict bound changes a pivot choice compared with exact minimum degree.
    offsets=[0,1,4,8,9,12,14,18,19,21,23,26,30,33,34,36,41,44,47,50]
    indices=[0,1,12,17,0,2,7,15,3,4,11,15,5,9,0,6,7,14,7,8,9,9,14,2,10,14,6,11,16,18,8,12,17,13,3,14,3,11,12,15,18,4,11,16,8,14,17,4,17,18]
    a=csr(19,19,offsets,indices,[1]*len(indices))
    cases.append(dict(name='strict bound changes pivot',op='amd',a=a,b=[0]*19,expected=amd(a),invalid=False))
    from cholesky_reference import fixtures as cholesky_cases
    for c in cholesky_cases():
        if c['op']=='chol_rcm_total':cases.append(dict(c,op='chol_amd_total',name=c['name'].replace('rcm','amd')))
    return cases
