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

    def approximate_minimum_degree(self):
        """Deterministic quotient-graph AMD; zero-based new-to-old permutation."""
        if self.rows != self.cols:
            raise ValueError('AMD requires square matrix')
        n = self.rows
        direct = [set() for _ in range(n)]
        elements = [set() for _ in range(n)]
        for i in range(n):
            for j in self._ci[self._rp[i]:self._rp[i+1]]:
                if i != j:
                    direct[i].add(j)
                    direct[j].add(i)
        active = set(range(n))
        degree = [len(g) for g in direct]
        order = []
        while active:
            pivot = min(active, key=lambda i: (min(degree[i], len(active)-1), i))
            neighbors = direct[pivot].copy()
            absorbed = [e for e in range(n) if pivot in elements[e]]
            for e in absorbed:
                neighbors.update(elements[e])
            neighbors.discard(pivot)
            active.remove(pivot)
            order.append(pivot)
            for e in absorbed:
                elements[e].clear()
            direct[pivot].clear()
            for i in neighbors:
                direct[i].discard(pivot)
                direct[i].difference_update(neighbors)
            elements[pivot] = neighbors
            for i in neighbors:
                bound = len(neighbors)-1 + len(direct[i])
                for e in range(n):
                    if e != pivot and i in elements[e]:
                        bound += len(elements[e] - neighbors)
                degree[i] = min(len(active)-1, bound)
        return order

    def reverse_cuthill_mckee(self):
        """New-to-old ordering of the undirected stored pattern, including zeros."""
        if self.rows!=self.cols:raise ValueError('RCM requires square matrix')
        graph=[set() for _ in range(self.rows)]
        for i in range(self.rows):
            for j in self._ci[self._rp[i]:self._rp[i+1]]:
                if i!=j:graph[i].add(j);graph[j].add(i)
        key=lambda i:(len(graph[i]),i)
        seen=set();order=[]
        for start in sorted(range(self.rows),key=key):
            if start in seen:continue
            queue=[start];seen.add(start)
            for i in queue:
                for j in sorted(graph[i]-seen,key=key):seen.add(j);queue.append(j)
            order.extend(queue)
        return order[::-1]

    @staticmethod
    def permute_vector(order,x,inverse=False):
        p=list(order);x=list(x);n=len(x)
        if len(p)!=n or any(type(i) is not int or i<0 or i>=n for i in p) or len(set(p))!=n or any(not math.isfinite(v) for v in x):raise ValueError('invalid permutation or vector')
        if not inverse:return [x[i] for i in p]
        out=[0.]*n
        for i,j in enumerate(p):out[j]=x[i]
        return out

    def permute_symmetric(self,order):
        p=list(order);n=self.rows
        if self.cols!=n:raise ValueError('permutation requires square matrix')
        inv=self.permute_vector(p,list(range(n)),True)
        rp,ci,v=[0],[],[]
        for i in p:
            for j,value in sorted((inv[self._ci[k]],self._v[k]) for k in range(self._rp[i],self._rp[i+1])):ci.append(j);v.append(value)
            rp.append(len(v))
        return CSRMatrix(n,n,rp,ci,v)

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

    def conjugate_gradient(self,b,rtol=1e-10,atol=0.0,max_iterations=1000,jacobi=False,capture=False,preconditioner=None):
        b=list(b);n=self.rows
        if self.cols!=n or len(b)!=n or any(not math.isfinite(v) for v in b):raise ValueError('CG requires square matrix and finite matching vector')
        if not math.isfinite(rtol) or not 0<=rtol<1 or not math.isfinite(atol) or atol<0 or type(max_iterations) is not int or not 0<=max_iterations<=100000:
            raise ValueError('invalid CG options')
        if preconditioner is not None:
            from cholesky import IC0
            from multigrid import GeometricMultigrid
            if not isinstance(preconditioner,(IC0,GeometricMultigrid)) or preconditioner.size!=n or jacobi:raise ValueError('invalid or conflicting CG preconditioner')
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
        def apply(r):return preconditioner.apply(r) if preconditioner is not None else [r[i]/diagonal[i] for i in range(n)]
        try:z=apply(r)
        except ArithmeticError:return result('nonfinite')
        p=z.copy();rho=dot(r,z)
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
            try:z=apply(r)
            except ArithmeticError:return result('nonfinite')
            next_rho=dot(r,z)
            if not math.isfinite(next_rho):return result('nonfinite')
            if next_rho<=0:return result('breakdown')
            beta=next_rho/rho;p=[z[i]+beta*p[i] for i in range(n)];rho=next_rho
        return result('iteration_limit')

    def gmres(self,b,restart=30,rtol=1e-10,atol=0.0,max_iterations=1000,jacobi=False,capture=False,preconditioner=None):
        from gmres import gmres
        return gmres(self,b,restart,rtol,atol,max_iterations,jacobi,capture,preconditioner)
