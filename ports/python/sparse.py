"""Canonical CSR storage and preconditioned conjugate gradient for real SPD systems."""
import math
from typing import NamedTuple


class CGResult(NamedTuple):
    x: list
    converged: bool
    iterations: int
    reason: str
    residuals: list
    iterates: list


class CSRMatrix:
    def __init__(self, rows, cols, row_offsets, column_indices, values):
        if type(rows) is not int or type(cols) is not int or rows < 0 or cols < 0:
            raise ValueError('invalid CSR dimensions')
        rp,ci,v=list(row_offsets),list(column_indices),[float(x) for x in values]
        if len(rp)!=rows+1 or len(ci)!=len(v) or rp[0]!=0 or rp[-1]!=len(v):
            raise ValueError('invalid CSR array lengths or endpoints')
        if any(type(x) is not int or x<0 or x>len(v) for x in rp) or any(not math.isfinite(x) for x in v):
            raise ValueError('invalid CSR offsets or values')
        for i in range(rows):
            if rp[i]>rp[i+1]:raise ValueError('CSR offsets must be monotone')
            previous=-1
            for p in range(rp[i],rp[i+1]):
                if type(ci[p]) is not int or not previous<ci[p]<cols:raise ValueError('CSR columns must be in range, sorted and unique')
                previous=ci[p]
        self._rows,self._cols,self._rp,self._ci,self._v=rows,cols,rp,ci,v

    rows=property(lambda self:self._rows)
    cols=property(lambda self:self._cols)
    nnz=property(lambda self:len(self._v))
    row_offsets=property(lambda self:self._rp.copy())
    column_indices=property(lambda self:self._ci.copy())
    values=property(lambda self:self._v.copy())

    @classmethod
    def from_dense(cls,a):
        rp,ci,v=[0],[],[]
        for i in range(a.rows):
            for j in range(a.cols):
                if a[i,j]!=0:ci.append(j);v.append(a[i,j])
            rp.append(len(v))
        return cls(a.rows,a.cols,rp,ci,v)

    def matvec(self,x):
        x=list(x)
        if len(x)!=self.cols or any(not math.isfinite(v) for v in x):raise ValueError('invalid vector')
        out=[sum(self._v[p]*x[self._ci[p]] for p in range(self._rp[i],self._rp[i+1])) for i in range(self.rows)]
        if any(not math.isfinite(v) for v in out):raise ArithmeticError('sparse multiplication outside float64 range')
        return out

    def conjugate_gradient(self,b,rtol=1e-10,atol=0.0,max_iterations=1000,jacobi=False,capture=False):
        b=list(b);n=self.rows
        if self.cols!=n or len(b)!=n or any(not math.isfinite(v) for v in b):raise ValueError('CG requires square matrix and finite matching vector')
        if not math.isfinite(rtol) or not 0<=rtol<1 or not math.isfinite(atol) or atol<0 or type(max_iterations) is not int or not 0<=max_iterations<=100000:
            raise ValueError('invalid CG options')
        # Validate symmetry without creating a dense copy. Explicit zero entries are allowed.
        from bisect import bisect_left
        diagonal=[1.0]*n
        for i in range(n):
            for p in range(self._rp[i],self._rp[i+1]):
                j=self._ci[p];q=bisect_left(self._ci,i,self._rp[j],self._rp[j+1])
                other=self._v[q] if q<self._rp[j+1] and self._ci[q]==i else 0
                if self._v[p]!=other:raise ValueError('CG requires exact symmetry')
            if jacobi:
                q=bisect_left(self._ci,i,self._rp[i],self._rp[i+1])
                if q==self._rp[i+1] or self._ci[q]!=i or self._v[q]<=0:raise ValueError('Jacobi requires a positive diagonal')
                diagonal[i]=self._v[q]
        def norm(v):
            s=0.0
            for a in v:s=math.hypot(s,a)
            return s
        def dot(a,b):return sum(x*y for x,y in zip(a,b))
        x=[0.0]*n;r=b.copy();history=[norm(r)];frames=[x.copy()] if capture else []
        threshold=max(atol,rtol*history[0])
        def result(reason):return CGResult(x.copy(),reason=='converged',len(history)-1,reason,history.copy(),frames.copy())
        if not math.isfinite(history[0]):return result('nonfinite')
        if history[0]<=threshold:return result('converged')
        z=[r[i]/diagonal[i] for i in range(n)];p=z.copy();rho=dot(r,z)
        for _ in range(max_iterations):
            try:q=self.matvec(p)
            except (ArithmeticError,ValueError):return result('nonfinite')
            curvature=dot(p,q)
            if not math.isfinite(rho) or not math.isfinite(curvature):return result('nonfinite')
            if rho<=0 or curvature<=0:return result('breakdown')
            alpha=rho/curvature
            candidate=[x[i]+alpha*p[i] for i in range(n)]
            try:ax=self.matvec(candidate)
            except (ArithmeticError,ValueError):return result('nonfinite')
            residual=[b[i]-ax[i] for i in range(n)];length=norm(residual)
            if not math.isfinite(length):return result('nonfinite')
            x,r=candidate,residual;history.append(length)
            if capture:frames.append(x.copy())
            if length<=threshold:return result('converged')
            z=[r[i]/diagonal[i] for i in range(n)];next_rho=dot(r,z)
            if not math.isfinite(next_rho):return result('nonfinite')
            if next_rho<=0:return result('breakdown')
            beta=next_rho/rho;p=[z[i]+beta*p[i] for i in range(n)];rho=next_rho
        return result('iteration_limit')
