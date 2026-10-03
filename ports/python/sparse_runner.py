import sys,math
from time import perf_counter_ns
from sparse import CSRMatrix
from matrix import Matrix

def run(args):
    if len(args)!=10:raise ValueError('expected sparse OP ROWS COLS NNZ ITERATIONS RTOL ATOL LIMIT JACOBI CAPTURE')
    op=args[0];rows,cols,nnz,iters=map(int,args[1:5]);rtol,atol=map(float,args[5:7]);limit,jacobi,capture=map(int,args[7:])
    if op not in ('spmv','dense','cg') or min(rows,cols,nnz,iters)<0 or jacobi not in (0,1) or capture not in (0,1):raise ValueError('invalid sparse protocol')
    raw=sys.stdin.read().split();count=cols if op!='cg' else rows
    if len(raw)!=rows+1+2*nnz+count:raise ValueError('incorrect sparse input count')
    rp=list(map(int,raw[:rows+1]));ci=list(map(int,raw[rows+1:rows+1+nnz]));v=list(map(float,raw[rows+1+nnz:rows+1+2*nnz]));b=list(map(float,raw[rows+1+2*nnz:]))
    a=CSRMatrix(rows,cols,rp,ci,v)
    dense=right=None
    if op=='dense':
        values=[0.]*(rows*cols)
        for i in range(rows):
            for p in range(rp[i],rp[i+1]):values[i*cols+ci[p]]=v[p]
        dense,right=Matrix(rows,cols,values),Matrix(cols,1,b)
    def compute():
        if op=='dense':return dense.multiply(right).values
        if op=='spmv':return a.matvec(b)
        r=a.conjugate_gradient(b,rtol,atol,limit,bool(jacobi),bool(capture))
        if iters:
            if not r.converged:raise ValueError('benchmark CG did not converge: '+r.reason)
            return r.x
        return [['converged','iteration_limit','breakdown','nonfinite'].index(r.reason),r.iterations,len(r.residuals)]+r.x+r.residuals+[x for frame in r.iterates for x in frame]
    if not iters:
        out=compute();return dict(rows=1,cols=len(out),values=out)
    for _ in range(3):compute()
    checksum=0.;start=perf_counter_ns()
    for _ in range(iters):checksum+=sum(compute())
    elapsed=perf_counter_ns()-start
    if not math.isfinite(checksum):raise ArithmeticError('nonfinite benchmark checksum')
    return dict(elapsed_ns=elapsed,iterations=iters,checksum=checksum)
