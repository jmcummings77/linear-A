"""Owned ILU(0): no fill, pivoting, reordering or diagonal shifts."""
import math
from bisect import bisect_left
class ILU0:
    def __init__(self,a):
        if a.rows!=a.cols:raise ValueError('ILU0 requires square matrix')
        self.size=a.rows;self._rp=list(a._rp);self._ci=list(a._ci);self._v=list(a._v);self._d=[]
        rp,ci,v=self._rp,self._ci,self._v
        for i in range(self.size):
            p=bisect_left(ci,i,rp[i],rp[i+1])
            if p==rp[i+1] or ci[p]!=i:raise ValueError('ILU0 requires stored diagonal')
            self._d.append(p)
        for i in range(self.size):
            for p in range(rp[i],self._d[i]):
                j=ci[p];v[p]/=v[self._d[j]]
                if not math.isfinite(v[p]):raise ArithmeticError('nonfinite ILU0 factor')
                for q in range(self._d[j]+1,rp[j+1]):
                    k=bisect_left(ci,ci[q],p+1,rp[i+1])
                    if k<rp[i+1] and ci[k]==ci[q]:
                        v[k]-=v[p]*v[q]
                        if not math.isfinite(v[k]):raise ArithmeticError('nonfinite ILU0 factor')
            if v[self._d[i]]==0:raise ArithmeticError('zero ILU0 pivot')
    @property
    def nnz(self):return len(self._v)
    def apply(self,b):
        x=list(b);rp,ci,v,d=self._rp,self._ci,self._v,self._d
        if len(x)!=self.size or any(not math.isfinite(z) for z in x):raise ValueError('invalid ILU0 vector')
        for i in range(self.size):
            for p in range(rp[i],d[i]):x[i]-=v[p]*x[ci[p]]
            if not math.isfinite(x[i]):raise ArithmeticError('nonfinite ILU0 solve')
        for i in range(self.size-1,-1,-1):
            for p in range(d[i]+1,rp[i+1]):x[i]-=v[p]*x[ci[p]]
            x[i]/=v[d[i]]
            if not math.isfinite(x[i]):raise ArithmeticError('nonfinite ILU0 solve')
        return x
