/** Independent browser gate: analytic solutions and true residuals for every frame. */
export function checkSparse(api,cases) {
  const norm=v=>v.reduce((s,x)=>Math.hypot(s,x),0);
  const assert=(condition,message)=>{if(!condition)throw new Error(message);};
  for(const test of cases) {
    let a,output,error;
    try {
      const m=test.a;a=new api.CSRMatrix(m.rows,m.cols,m.offsets,m.indices,m.values);
      if(test.op==='cg'||test.op==='gmres') {const o=test.options,options={restart:o.restart,rtol:o.rtol,atol:o.atol,maxIterations:o.limit,jacobi:!!o.jacobi,capture:!!o.capture};output=test.op==='gmres'?a.gmres(test.b,options):a.conjugateGradient(test.b,options);}
      else if(test.op==='spmv') output=Array.from(a.matvec(test.b));
      else {
        const v=Array(m.rows*m.cols).fill(0);for(let i=0;i<m.rows;i++)for(let p=m.offsets[i];p<m.offsets[i+1];p++)v[i*m.cols+m.indices[p]]=m.values[p];
        const dense=new api.Matrix(m.rows,m.cols,v);let b,c;
        try{b=new api.Matrix(m.cols,1,test.b);c=dense.multiply(b);output=Array.from(c.toArray?c.toArray():c.values);}
        finally{dense.dispose?.();b?.dispose?.();c?.dispose?.();}
      }
    } catch(e){error=e;} finally{a?.dispose?.();}
    if(test.invalid){assert(error,`${test.name}: accepted invalid input`);continue;}
    if(error)throw new Error(`${test.name}: ${error.message}`);
    const close=(x,y)=>x.length===y.length&&x.every((v,i)=>Number.isFinite(v)&&Math.abs(v-y[i])<=1e-7*Math.max(1,Math.abs(y[i])));
    if(test.op!=='cg'&&test.op!=='gmres'){assert(close(output,test.expected),test.name);continue;}
    const r=output,n=test.a.rows;
    assert(r.reason===test.reason&&r.converged===(r.reason==='converged'),`${test.name}: stop reason`);
    assert(Number.isInteger(r.iterations)&&r.iterations>=0&&r.iterations<=test.options.limit&&r.residuals.length===r.iterations+1,`${test.name}: history length`);
    assert(r.x.length===n&&r.x.every(Number.isFinite)&&r.residuals.every(Number.isFinite),`${test.name}: finite result`);
    if(test.expected)assert(close(r.x,test.expected),`${test.name}: solution`);
    const residual=x=>norm(test.b.map((b,i)=>{let sum=0;for(let p=test.a.offsets[i];p<test.a.offsets[i+1];p++)sum+=test.a.values[p]*x[test.a.indices[p]];return b-sum;}));
    if(r.converged)assert(residual(r.x)<=Math.max(test.options.atol,test.options.rtol*norm(test.b))*(1+1e-6)+1e-14,`${test.name}: convergence`);
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
