"""Restarted, right-preconditioned GMRES with two-pass modified Gram–Schmidt."""
import math
from typing import NamedTuple


class GMRESResult(NamedTuple):
    x: list
    converged: bool
    iterations: int
    reason: str
    residuals: list
    estimated_residuals: list
    iterates: list
    restarts: list


def gmres(a, b, restart=30, rtol=1e-10, atol=0.0, max_iterations=1000, jacobi=False, capture=False, preconditioner=None):
    b=list(b);n=a.rows
    if a.cols!=n or len(b)!=n or any(not math.isfinite(x) for x in b):raise ValueError('GMRES requires square matrix and finite matching vector')
    if type(restart) is not int or not 1<=restart<=1024 or type(max_iterations) is not int or not 0<=max_iterations<=100000 or not math.isfinite(rtol) or not 0<=rtol<1 or not math.isfinite(atol) or atol<0 or type(jacobi) is not bool or type(capture) is not bool:raise ValueError('invalid GMRES options')
    if preconditioner is not None and (jacobi or preconditioner.size!=n):raise ValueError('incompatible preconditioner')
    apply=lambda v:preconditioner.apply(v) if preconditioner is not None else [v[i]/diagonal[i] for i in range(n)]
    diagonal=[1.0]*n
    if jacobi:
        for i in range(n):
            for p in range(a._rp[i],a._rp[i+1]):
                if a._ci[p]==i:diagonal[i]=a._v[p];break
            else:raise ValueError('Jacobi requires a nonzero diagonal')
            if diagonal[i]==0:raise ValueError('Jacobi requires a nonzero diagonal')
    def norm(v):
        result=0.0
        for value in v:result=math.hypot(result,value)
        return result
    def finite(v):return all(math.isfinite(x) for x in v)
    x=[0.0]*n;r=b.copy();history=[norm(r)];estimates=history.copy();frames=[x.copy()] if capture else [];restarts=[]
    def result(reason):return GMRESResult(x.copy(),reason=='converged',len(history)-1,reason,history.copy(),estimates.copy(),frames.copy(),restarts.copy())
    threshold=max(atol,rtol*history[0]);m=min(restart,n,max_iterations)
    if not math.isfinite(history[0]):return result('nonfinite')
    if history[0]<=threshold:return result('converged')
    while len(history)-1<max_iterations:
        if len(history)>1:restarts.append(len(history)-1)
        base=x.copy();beta=norm(r);basis=[[v/beta for v in r]]
        h=[[0.0]*m for _ in range(m+1)];cs=[0.0]*m;sn=[0.0]*m;g=[beta]+[0.0]*m
        for j in range(min(m,max_iterations-(len(history)-1))):
            try:w=a.matvec(apply(basis[j]))
            except (ValueError,ArithmeticError):return result('nonfinite')
            original=norm(w)
            # A second MGS pass limits loss of orthogonality on difficult systems.
            for _ in range(2):
                for k in range(j+1):
                    dot=sum(basis[k][i]*w[i] for i in range(n));h[k][j]+=dot
                    w=[w[i]-dot*basis[k][i] for i in range(n)]
            tail=norm(w)
            if not math.isfinite(original) or not math.isfinite(tail) or any(not math.isfinite(h[k][j]) for k in range(j+1)):return result('nonfinite')
            happy=tail<=8*2.220446049250313e-16*original
            h[j+1][j]=0.0 if happy else tail
            if not happy:basis.append([v/tail for v in w])
            for k in range(j):
                top=cs[k]*h[k][j]+sn[k]*h[k+1][j]
                h[k+1][j]=-sn[k]*h[k][j]+cs[k]*h[k+1][j];h[k][j]=top
            pivot=math.hypot(h[j][j],h[j+1][j])
            if not math.isfinite(pivot):return result('nonfinite')
            if pivot==0:return result('breakdown')
            cs[j]=h[j][j]/pivot;sn[j]=h[j+1][j]/pivot;h[j][j]=pivot;h[j+1][j]=0.0
            g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];y=g[:j+1]
            for k in range(j,-1,-1):
                if h[k][k]==0:return result('breakdown')
                y[k]=(y[k]-sum(h[k][q]*y[q] for q in range(k+1,j+1)))/h[k][k]
            try:correction=apply([sum(basis[k][i]*y[k] for k in range(j+1)) for i in range(n)])
            except (ValueError,ArithmeticError):return result('nonfinite')
            candidate=[base[i]+correction[i] for i in range(n)]
            if not finite(y) or not finite(candidate) or not math.isfinite(g[j+1]):return result('nonfinite')
            try:ax=a.matvec(candidate)
            except (ValueError,ArithmeticError):return result('nonfinite')
            residual=[b[i]-ax[i] for i in range(n)];length=norm(residual)
            if not math.isfinite(length):return result('nonfinite')
            x,r=candidate,residual;history.append(length);estimates.append(abs(g[j+1]))
            if capture:frames.append(x.copy())
            if length<=threshold:return result('converged')
            if happy:return result('breakdown')
        if x==base:return result('stagnation')
    return result('iteration_limit')
