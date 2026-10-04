"""Symmetric geometric V-cycle for the unit five-point Dirichlet Laplacian."""
import math
from sparse import CSRMatrix

class GeometricMultigrid:
    """Nested interior grids of width 2**k-1; fixed two-sweep weighted Jacobi."""
    def __init__(self,width):
        if type(width) is not int or not 1<=width<=255 or (width+1)&width:raise ValueError('grid width must be 2**k-1 in 1..255')
        self._width=width
    width=property(lambda self:self._width)
    size=property(lambda self:self._width**2)
    levels=property(lambda self:(self._width+1).bit_length()-1)
    @staticmethod
    def _mv(w,x):
        out=[]
        for i,v in enumerate(x):
            y,j=divmod(i,w);s=4*v
            if j:s-=x[i-1]
            if j+1<w:s-=x[i+1]
            if y:s-=x[i-w]
            if y+1<w:s-=x[i+w]
            out.append(s)
        return out
    @property
    def matrix(self):
        w=self.width;rp=[0];ci=[];v=[]
        for i in range(self.size):
            y,j=divmod(i,w)
            for k in sorted([i]+([i-1] if j else [])+([i+1] if j+1<w else [])+([i-w] if y else [])+([i+w] if y+1<w else [])):ci.append(k);v.append(4. if k==i else -1.)
            rp.append(len(v))
        return CSRMatrix(self.size,self.size,rp,ci,v)
    def apply(self,b):return self._run(b,False)['x']
    def trace(self,b):return self._run(b,True)
    def _run(self,b,capture):
        b=list(b)
        if len(b)!=self.size or any(not math.isfinite(v) for v in b):raise ValueError('invalid multigrid RHS')
        frames=[]
        def record(level,w,phase,x,b):
            if any(not math.isfinite(v) for v in x):raise ArithmeticError('nonfinite multigrid cycle')
            if capture:frames.append(dict(level=level,width=w,phase=phase,x=x.copy(),b=b.copy(),residual=[v-a for v,a in zip(b,self._mv(w,x))]))
        def cycle(w,b,level):
            x=[0.]*len(b);record(level,w,'enter',x,b)
            if w==1:
                x=[b[0]/4];record(level,w,'coarse_solve',x,b);return x
            def smooth(x):
                for _ in range(2):x=[a+(v-q)/6 for a,v,q in zip(x,b,self._mv(w,x))]
                return x
            x=smooth(x);record(level,w,'pre_smooth',x,b);r=[v-a for v,a in zip(b,self._mv(w,x))];c=w//2;bc=[0.]*(c*c)
            for y in range(c):
                for j in range(c):
                    for dy in (-1,0,1):
                        for dx in (-1,0,1):bc[y*c+j]+=(1 if dy==0 else .5)*(1 if dx==0 else .5)*r[(2*y+1+dy)*w+2*j+1+dx]
            ec=cycle(c,bc,level+1)
            for y in range(c):
                for j in range(c):
                    for dy in (-1,0,1):
                        for dx in (-1,0,1):x[(2*y+1+dy)*w+2*j+1+dx]+=(1 if dy==0 else .5)*(1 if dx==0 else .5)*ec[y*c+j]
            record(level,w,'correct',x,b);x=smooth(x);record(level,w,'post_smooth',x,b);return x
        x=cycle(self.width,b,0)
        return dict(x=x,frames=frames)
