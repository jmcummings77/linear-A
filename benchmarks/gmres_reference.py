"""Independent GMRES fixtures and upwind advection–diffusion assembly."""
import math
from sparse_reference import csr,dense_csr,multiply,norm


def transport(size,contrast=0,diffusivity=1.0,speed=4.0,angle=30.0):
    a=diffusion_scaled(size,contrast,diffusivity)
    vx=speed*math.cos(math.radians(angle))/(size+1)
    vy=speed*math.sin(math.radians(angle))/(size+1)
    rp,ci,values=[0],[],[]
    for i in range(size*size):
        y,x=divmod(i,size);entries={a['indices'][p]:a['values'][p] for p in range(a['offsets'][i],a['offsets'][i+1])}
        for component,dy,dx in [(vx,0,-1 if vx>=0 else 1),(vy,-1 if vy>=0 else 1,0)]:
            weight=abs(component);entries[i]+=weight;yy,xx=y+dy,x+dx
            if 0<=yy<size and 0<=xx<size:
                j=yy*size+xx;entries[j]=entries.get(j,0)-weight
        for j,value in sorted(entries.items()):ci.append(j);values.append(value)
        rp.append(len(values))
    return csr(size*size,size*size,rp,ci,values)


def diffusion_scaled(size,contrast,diffusivity):
    from sparse_reference import diffusion
    a=diffusion(size,contrast);a['values']=[v*diffusivity for v in a['values']];return a


def dense_solve(a,b):
    """Small partial-pivot elimination oracle, separate from every sparse algorithm."""
    n=a['rows'];rows=[[0.0]*n+[b[i]] for i in range(n)]
    for i in range(n):
        for p in range(a['offsets'][i],a['offsets'][i+1]):rows[i][a['indices'][p]]=a['values'][p]
    for k in range(n):
        pivot=max(range(k,n),key=lambda i:abs(rows[i][k]));rows[k],rows[pivot]=rows[pivot],rows[k]
        if rows[k][k]==0:raise ValueError('singular oracle matrix')
        for i in range(k+1,n):
            factor=rows[i][k]/rows[k][k]
            for j in range(k+1,n+1):rows[i][j]-=factor*rows[k][j]
    x=[0.0]*n
    for i in range(n-1,-1,-1):x[i]=(rows[i][n]-sum(rows[i][j]*x[j] for j in range(i+1,n)))/rows[i][i]
    return x


def fixtures():
    cases=[]
    def add(name,a,b,expected=None,reason='converged',invalid=False,**options):
        cases.append(dict(name='GMRES '+name,a=a,b=b,op='gmres',expected=expected,reason=reason,invalid=invalid,
                          options={**dict(restart=3,rtol=1e-10,atol=0,limit=100,jacobi=0,capture=1),**options}))
    a=dense_csr([[3,1],[0,2]])
    add('nonsymmetric analytic',a,[5,4],expected=[1,2])
    add('right Jacobi',a,[5,4],expected=[1,2],jacobi=1)
    add('restart one',a,[5,4],expected=[1,2],restart=1)
    add('partial last cycle',dense_csr([[4,1,0],[0,3,1],[1,0,2]]),[1,2,3],reason='iteration_limit',restart=2,limit=3)
    add('no iterations',a,[5,4],reason='iteration_limit',limit=0)
    add('absolute stop',a,[1e-12,0],expected=[0,0],atol=1e-11)
    add('zero RHS',a,[0,0],expected=[0,0])
    add('empty',csr(0,0,[0],[],[]),[],expected=[])
    add('no capture',a,[5,4],expected=[1,2],capture=0)
    add('large restart capped by dimension',a,[5,4],expected=[1,2],restart=1024)
    add('happy breakdown identity',dense_csr([[1,0],[0,1]]),[2,-3],expected=[2,-3])
    add('negative Jacobi diagonal',dense_csr([[-2,1],[0,3]]),[0,6],expected=[1,2],jacobi=1)
    rotation=dense_csr([[0,-1],[1,0]])
    add('rotation indefinite',rotation,[-2,1],expected=[1,2],restart=2)
    add('restart stagnation',rotation,[1,0],reason='stagnation',restart=1)
    add('zero operator breakdown',dense_csr([[0,0],[0,0]]),[1,2],reason='breakdown')
    add('singular inconsistent',dense_csr([[1,0],[0,0]]),[0,1],reason='breakdown')
    add('false happy breakdown',dense_csr([[1,0],[0,0]]),[1,1],reason='breakdown')
    add('tiny RHS',dense_csr([[1]]),[1e-200],expected=[1e-200])
    add('nonfinite product',dense_csr([[1.7e308,1.7e308],[1.7e308,1.7e308]]),[1,1],reason='nonfinite')
    add('preconditioner overflow',dense_csr([[1e-320]]),[1],jacobi=1,reason='nonfinite')
    add('missing Jacobi diagonal',rotation,[1,2],jacobi=1,invalid=True)
    add('zero Jacobi diagonal',csr(1,1,[0,1],[0],[0]),[1],jacobi=1,invalid=True)
    add('rectangular',csr(1,2,[0,1],[0],[1]),[1],invalid=True)
    for key,value in [('restart',0),('restart',1025),('rtol',-1),('rtol',1),('atol',-1),('limit',100001)]:add('invalid '+key+str(value),a,[5,4],invalid=True,**{key:value})
    for restart in [2,5]:
        for jacobi in [0,1]:
            a=transport(3,1,.2,3,35);truth=[math.sin(i+1) for i in range(9)];b=multiply(a,truth)
            reference=dense_solve(a,b)
            assert max(abs(x-y) for x,y in zip(reference,truth))<1e-12
            add('transport restart %d Jacobi %d'%(restart,jacobi),a,b,expected=reference,restart=restart,jacobi=jacobi,limit=300)
    return cases


def unpack(actual,n,capture):
    v=actual['values'];reason,iterations,length=map(int,v[:3]);x=v[3:3+n];history=v[3+n:3+n+length];offset=3+n+length
    estimates=v[offset:offset+length];offset+=length
    count=v[offset];assert int(count)==count and count>=0;count=int(count);offset+=1
    restarts=v[offset:offset+count];flat=v[offset+count:]
    assert len(flat)==(n*length if capture else 0)
    return dict(reason=reason,iterations=iterations,x=x,residuals=history,estimated_residuals=estimates,restarts=restarts,iterates=[flat[k*n:(k+1)*n] for k in range(length)] if capture else [])


def check_result(actual,case):
    assert set(actual)=={'rows','cols','values'} and actual['rows']==1 and actual['cols']==len(actual['values'])
    assert all(type(v) in (float,int) and math.isfinite(v) for v in actual['values'])
    assert all(int(v)==v for v in actual['values'][:3])
    n=case['a']['rows'];o=case['options'];r=unpack(actual,n,o['capture']);length=r['iterations']+1
    assert 0<=r['reason']<5 and ['converged','iteration_limit','breakdown','nonfinite','stagnation'][r['reason']]==case['reason']
    assert actual['values'][2]==length and 0<=r['iterations']<=o['limit']
    assert len(r['x'])==n and len(r['residuals'])==len(r['estimated_residuals'])==length
    assert all(v>=0 for v in r['residuals']+r['estimated_residuals'])
    assert r['residuals'][0]==r['estimated_residuals'][0]
    m=min(o['restart'],n,o['limit']);expected_starts=list(range(m,r['iterations'],m)) if m else []
    # A cycle may fail before accepting its first update: that start is retained.
    assert r['restarts']==expected_starts or (r['reason'] in (2,3) and r['restarts']==expected_starts+[r['iterations']] and r['iterations']>0 and r['iterations']%m==0)
    if case['expected'] is not None:assert all(abs(x-y)<=1e-7*max(abs(y),1e-200) for x,y in zip(r['x'],case['expected']))
    residual=lambda x:norm([b-v for b,v in zip(case['b'],multiply(case['a'],x))])
    if r['reason']==0:assert residual(r['x'])<=max(o['atol'],o['rtol']*norm(case['b']))*(1+1e-6)+1e-300
    for k,frame in enumerate(r['iterates']):
        value=residual(frame);assert abs(value-r['residuals'][k])<=1e-11*max(value,r['residuals'][0],1e-300),'incorrect true residual'
    if r['iterates']:assert r['x']==r['iterates'][-1]
    if case['name']=='GMRES nonsymmetric analytic':
        expected=[635/425,508/425]
        assert all(abs(x-y)<1e-12 for x,y in zip(r['iterates'][1],expected))
        assert abs(r['estimated_residuals'][1]-residual(expected))<1e-12
    return r
