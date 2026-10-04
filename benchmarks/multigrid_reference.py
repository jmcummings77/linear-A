"""Independent dense transfer-matrix oracle for the structured V-cycle.

Small fixtures use exact rational arithmetic; benchmark grids use float64.
Restriction is the transpose of bilinear interpolation, and every level uses
its rediscretized, unscaled five-point operator (not a Galerkin product).
"""
from fractions import Fraction
import math
from sparse_reference import csr, multiply


def poisson(width):
    if type(width) is not int or not 1 <= width <= 255 or width & (width + 1):
        raise ValueError('width must be 2^k-1 in 1..255')
    n=width*width;rp=[0];ci=[];v=[]
    for i in range(n):
        for j in sorted({i} | ({i-width} if i>=width else set()) | ({i+width} if i+width<n else set()) | ({i-1} if i%width else set()) | ({i+1} if i%width+1<width else set())):
            ci.append(j);v.append(4. if i==j else -1.)
        rp.append(len(v))
    return csr(n,n,rp,ci,v)


def transfer(width):
    """Dense P from tensor products of one-dimensional hat functions."""
    c=width//2
    hats=[[max(Fraction(0),1-abs(Fraction(i+1,2)-(j+1))) for j in range(c)] for i in range(width)]
    return [[hats[i//width][j//c]*hats[i%width][j%c] for j in range(c*c)] for i in range(width*width)]


def reference(width,b):
    """Exact dense oracle, deliberately limited to small conformance grids."""
    if width>7:raise ValueError('dense rational oracle is for widths up to seven')
    a=poisson(width);n=len(b)
    if n!=a['rows']:raise ValueError('RHS shape')
    b=list(map(lambda v:Fraction(str(v)),b))
    dense=[[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for p in range(a['offsets'][i],a['offsets'][i+1]):dense[i][a['indices'][p]]=Fraction(a['values'][p])
    def mv(m,x):return [sum(v*y for v,y in zip(row,x)) for row in m]
    if width==1:return [b[0]/4]
    x=[Fraction(0)]*n
    def smooth(x):
        for _ in range(2):x=[v+(r-q)/6 for v,r,q in zip(x,b,mv(dense,x))]
        return x
    x=smooth(x);r=[v-q for v,q in zip(b,mv(dense,x))];p=transfer(width)
    coarse=mv(list(zip(*p)),r);correction=reference(width//2,coarse)
    x=[v+q for v,q in zip(x,mv(p,correction))]
    return smooth(x)


def truth(width):
    # Deterministic broad spectrum: prevents CG's early termination on a single eigenmode.
    return [math.sin((i+1)*1.731)+0.2*math.cos((i+3)*.413) for i in range(width*width)]


def fixtures():
    cases=[]
    def add(name,a,b,op,expected=None,invalid=False,**opts):
        cases.append(dict(name='Multigrid '+name,a=a,b=b,op=op,expected=expected,invalid=invalid,reason='converged',options=dict(rtol=1e-10,atol=0,limit=200,capture=int(op=='cg'),jacobi=opts.get('jacobi',0))))
    for w in (1,3,7):
        a=poisson(w);b=[float((i*7)%11-5) for i in range(w*w)]
        add('setup '+str(w),a,b,'mg_setup',[w*w,(w+1).bit_length()-1])
        add('matrix '+str(w),a,b,'mg_matrix',a['offsets']+a['indices']+a['values'])
        add('apply '+str(w),a,b,'mg_apply',list(map(float,reference(w,b))))
        for mode in (4,5):add('CG '+str(w)+' mode '+str(mode),a,multiply(a,truth(w)),'cg',truth(w),jacobi=mode)
    a=poisson(3)
    add('zero RHS',a,[0.]*9,'mg_apply',[0.]*9)
    for n in (0,2,4,16):
        a=csr(n,n,list(range(n+1)),list(range(n)),[4.]*n)
        add('invalid shape '+str(n),a,[1.]*n,'mg_setup',invalid=True)
    add('overflow',poisson(3),[1.7e308]*9,'mg_apply',invalid=True)
    return cases
