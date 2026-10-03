import sys,math
from time import perf_counter_ns
from sparse import CSRMatrix
from ilu import ILU0
from matrix import Matrix

def run(args):
    if len(args)!=(11 if args and args[0]=='gmres' else 10):raise ValueError('expected sparse OP ROWS COLS NNZ ITERATIONS RTOL ATOL LIMIT JACOBI CAPTURE')
    op=args[0];rows,cols,nnz,iters=map(int,args[1:5]);rtol,atol=map(float,args[5:7]);limit,jacobi,capture=map(int,args[7:10]);restart=int(args[10]) if op=='gmres' else 30
    if op not in ('spmv','dense','cg','gmres','ilu_setup','ilu_apply','rcm','permute','permutation_check','rcm_solve','ilu_solve') or min(rows,cols,nnz,iters)<0 or jacobi not in (range(4) if op=='gmres' else (0,1)) or capture not in (0,1):raise ValueError('invalid sparse protocol')
    raw=sys.stdin.read().split();count=rows if op in ('cg','gmres') else cols
    if len(raw)!=rows+1+2*nnz+count:raise ValueError('incorrect sparse input count')
    rp=list(map(int,raw[:rows+1]));ci=list(map(int,raw[rows+1:rows+1+nnz]));v=list(map(float,raw[rows+1+nnz:rows+1+2*nnz]));b=list(map(float,raw[rows+1+2*nnz:]))
    a=CSRMatrix(rows,cols,rp,ci,v)
    factor=ILU0(a) if op=='ilu_apply' or (op=='gmres' and jacobi==2) else None
    dense=right=None
    if op=='dense':
        values=[0.]*(rows*cols)
        for i in range(rows):
            for p in range(rp[i],rp[i+1]):values[i*cols+ci[p]]=v[p]
        dense,right=Matrix(rows,cols,values),Matrix(cols,1,b)
    def compute():
        if op in ('rcm_solve','ilu_solve'):
            p=a.reverse_cuthill_mckee() if op=='rcm_solve' else list(range(rows))
            q=a.permute_symmetric(p) if op=='rcm_solve' else a
            rhs=CSRMatrix.permute_vector(p,b) if op=='rcm_solve' else b
            result=q.gmres(rhs,20,rtol,atol,limit,False,False,ILU0(q))
            if not result.converged:raise ArithmeticError('ordering solve failed: '+result.reason)
            x=CSRMatrix.permute_vector(p,result.x,True) if op=='rcm_solve' else result.x
            return [result.iterations]+x
        if op=='rcm':return a.reverse_cuthill_mckee()
        if op in ('permute','permutation_check'):
            if any(not math.isfinite(x) or x!=int(x) for x in b):raise ValueError('invalid permutation')
            p=list(map(int,b));q=a.permute_symmetric(p)
            if op=='permute':return q.row_offsets+q.column_indices+q.values
            x=[i+1. for i in range(rows)];y=CSRMatrix.permute_vector(p,x)
            return y+CSRMatrix.permute_vector(p,y,True)+CSRMatrix.permute_vector(p,q.matvec(y),True)
        if op=='ilu_setup':
            f=ILU0(a);return [f.size,f.nnz]
        if op=='ilu_apply':return factor.apply(b)
        if op=='dense':return dense.multiply(right).values
        if op=='spmv':return a.matvec(b)
        r=a.gmres(b,restart,rtol,atol,limit,jacobi==1,bool(capture),ILU0(a) if jacobi==3 else factor) if op=='gmres' else a.conjugate_gradient(b,rtol,atol,limit,bool(jacobi),bool(capture))
        if iters:
            if not r.converged:raise ValueError('benchmark CG did not converge: '+r.reason)
            return r.x
        return [['converged','iteration_limit','breakdown','nonfinite','stagnation'].index(r.reason),r.iterations,len(r.residuals)]+r.x+r.residuals+((r.estimated_residuals+[len(r.restarts)]+r.restarts) if op=='gmres' else [])+[x for frame in r.iterates for x in frame]
    if not iters:
        out=compute();return dict(rows=1,cols=len(out),values=out)
    for _ in range(3):compute()
    checksum=0.;start=perf_counter_ns()
    for _ in range(iters):checksum+=sum(compute())
    elapsed=perf_counter_ns()-start
    if not math.isfinite(checksum):raise ArithmeticError('nonfinite benchmark checksum')
    return dict(elapsed_ns=elapsed,iterations=iters,checksum=checksum)
