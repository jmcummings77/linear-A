/** Independent browser gate: analytic solutions and true residuals for every frame. */
export function checkSparse(api,cases) {
  const norm=v=>v.reduce((s,x)=>Math.hypot(s,x),0);
  const assert=(condition,message)=>{if(!condition)throw new Error(message);};
  for(const test of cases) {
    let a,factor,output,error;
    try {
      const m=test.a;a=new api.CSRMatrix(m.rows,m.cols,m.offsets,m.indices,m.values);
      if(test.op.startsWith('chol_')){
        let plan,l,q=a;const reorder=test.op==='chol_rcm_total'||test.op==='chol_amd_total';
        try{const p=test.op==='chol_amd_total'?a.approximateMinimumDegree():reorder?a.reverseCuthillMcKee():undefined;if(reorder)q=a.permuteSymmetric(p);plan=new api.SparseCholeskySymbolic(q);
          if(test.op==='chol_symbolic')output=[...plan.rowOffsets,...plan.columnIndices,...plan.fillSteps];
          else{factor=plan.factorize(q);if(test.op==='chol_factor'){l=factor.lower;output=[...l.rowOffsets,...l.columnIndices,...l.values];}
          else{const x=factor.solve(reorder?api.CSRMatrix.permuteVector(p,test.b):test.b);output=reorder?api.CSRMatrix.permuteVector(p,x,true):x;
            const ax=a.matvec(output);assert(norm(test.b.map((v,i)=>v-ax[i]))<=1e-11*Math.max(norm(test.b),1e-300),'Cholesky original-system residual');}}
        }finally{l?.dispose?.();plan?.dispose?.();if(q!==a)q.dispose?.();}
      }
      else if(test.op==='rcm_solve'||test.op==='ilu_solve'){
        const reorder=test.op==='rcm_solve',p=reorder?a.reverseCuthillMcKee():[],q=reorder?a.permuteSymmetric(p):a;
        try{factor=new api.ILU0(q);const rhs=reorder?api.CSRMatrix.permuteVector(p,test.b):test.b,r=q.gmres(rhs,{restart:20,rtol:test.options.rtol,atol:test.options.atol,maxIterations:test.options.limit,preconditioner:factor});assert(r.converged,'ordering solve failed');output=reorder?api.CSRMatrix.permuteVector(p,r.x,true):r.x;
        const ax=a.matvec(output);assert(norm(test.b.map((v,i)=>v-ax[i]))<=Math.max(test.options.atol,test.options.rtol*norm(test.b))*(1+1e-5)+1e-300,'original-system residual');}
        finally{if(reorder)q.dispose?.();}
      }
      else if(test.op==='amd')output=a.approximateMinimumDegree();
      else if(test.op==='rcm')output=a.reverseCuthillMcKee();
      else if(test.op==='permute'||test.op==='permutation_check'){
        const p=test.b,q=a.permuteSymmetric(p);
        try{if(test.op==='permute')output=[...q.rowOffsets,...q.columnIndices,...q.values];
        else{const x=Array.from({length:m.rows},(_,i)=>i+1),y=api.CSRMatrix.permuteVector(p,x);output=[...y,...api.CSRMatrix.permuteVector(p,y,true),...api.CSRMatrix.permuteVector(p,q.matvec(y),true)];}}
        finally{q.dispose?.();}
      }
      else if(test.op==='ilu_setup'){factor=new api.ILU0(a);output=[factor.size,factor.nnz];}
      else if(test.op==='ilu_apply'){factor=new api.ILU0(a);output=Array.from(factor.apply(test.b));}
      else if(test.op==='cg'||test.op==='gmres') {const o=test.options,options={restart:o.restart,rtol:o.rtol,atol:o.atol,maxIterations:o.limit,jacobi:o.jacobi===1,capture:!!o.capture};if(test.op==='gmres'&&o.jacobi>=2){factor=new api.ILU0(a);options.preconditioner=factor;}output=test.op==='gmres'?a.gmres(test.b,options):a.conjugateGradient(test.b,options);}
      else if(test.op==='spmv') output=Array.from(a.matvec(test.b));
      else {
        const v=Array(m.rows*m.cols).fill(0);for(let i=0;i<m.rows;i++)for(let p=m.offsets[i];p<m.offsets[i+1];p++)v[i*m.cols+m.indices[p]]=m.values[p];
        const dense=new api.Matrix(m.rows,m.cols,v);let b,c;
        try{b=new api.Matrix(m.cols,1,test.b);c=dense.multiply(b);output=Array.from(c.toArray?c.toArray():c.values);}
        finally{dense.dispose?.();b?.dispose?.();c?.dispose?.();}
      }
    } catch(e){error=e;} finally{factor?.dispose?.();a?.dispose?.();}
    if(test.invalid){assert(error,`${test.name}: accepted invalid input`);continue;}
    if(error)throw new Error(`${test.name}: ${error.message}`);
    const close=(x,y)=>x.length===y.length&&x.every((v,i)=>Number.isFinite(v)&&Math.abs(v-y[i])<=(test.op.startsWith('chol_')?1e-10:1e-7)*Math.max(test.op.startsWith('chol_')?1e-300:test.op==='gmres'?1e-200:1,Math.abs(y[i])));
    if(test.op!=='cg'&&test.op!=='gmres'){assert(close(output,test.expected),test.name);continue;}
    const r=output,n=test.a.rows;
    assert(r.reason===test.reason&&r.converged===(r.reason==='converged'),`${test.name}: stop reason`);
    assert(Number.isInteger(r.iterations)&&r.iterations>=0&&r.iterations<=test.options.limit&&r.residuals.length===r.iterations+1,`${test.name}: history length`);
    assert(r.x.length===n&&r.x.every(Number.isFinite)&&r.residuals.every(Number.isFinite),`${test.name}: finite result`);
    if(test.expected)assert(close(r.x,test.expected),`${test.name}: solution`);
    const residual=x=>norm(test.b.map((b,i)=>{let sum=0;for(let p=test.a.offsets[i];p<test.a.offsets[i+1];p++)sum+=test.a.values[p]*x[test.a.indices[p]];return b-sum;}));
    if(r.converged)assert(residual(r.x)<=Math.max(test.options.atol,test.options.rtol*norm(test.b))*(1+1e-6)+(test.op==='gmres'?1e-300:1e-14),`${test.name}: convergence`);
    assert(r.iterates.length===(test.options.capture?r.residuals.length:0),`${test.name}: frames`);
    r.iterates.forEach((frame,i)=>assert(frame.length===n&&Math.abs(residual(frame)-r.residuals[i])<=1e-11*Math.max(r.residuals[0],r.residuals[i],1e-300),`${test.name}: true residual`));
    if(r.iterates.length)assert(close(r.x,r.iterates.at(-1)),`${test.name}: final frame`);
    if(test.op==='gmres'){
      assert(r.estimatedResiduals.length===r.residuals.length&&r.estimatedResiduals.every(v=>Number.isFinite(v)&&v>=0),`${test.name}: estimates`);
      assert(r.estimatedResiduals[0]===r.residuals[0],`${test.name}: initial estimate`);
      const m=Math.min(test.options.restart,n,test.options.limit);
      r.restarts.forEach((v,i)=>assert(Number.isInteger(v)&&v===(i+1)*m&&v<=r.iterations,`${test.name}: restart boundaries`));
      assert(r.restarts.length===Math.floor(Math.max(0,r.iterations-1)/Math.max(1,m))||(r.restarts.at(-1)===r.iterations&&['breakdown','nonfinite'].includes(r.reason)),`${test.name}: restart count`);
      if(test.name==='GMRES nonsymmetric analytic'){
        const expected=[635/425,508/425];assert(close(r.iterates[1],expected),`${test.name}: first projection`);
        assert(Math.abs(residual(expected)-r.estimatedResiduals[1])<1e-12,`${test.name}: projected residual`);
      }
    }

  }
  return cases.length;
}
