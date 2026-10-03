/** Small, real WASM solves. Residuals are recomputed independently in JavaScript. */
function entries(values, count, name) {
  if (!Array.isArray(values) || values.length !== count || values.some(v => typeof v !== 'number' || !Number.isFinite(v) || Math.abs(v) > 1e6))
    throw new RangeError(`${name} must contain ${count} finite entries in [−10⁶, 10⁶]`);
  return values.slice();
}
const infinity = values => Math.max(0, ...values.map(Math.abs));

export function diagnostics(a, b, x) {
  const anorm = Math.max(Math.abs(a[0])+Math.abs(a[1]), Math.abs(a[2])+Math.abs(a[3]));
  const bnorm = infinity(b), xnorm = infinity(x), scale = Math.max(1, xnorm);
  // Divide by a common solution scale before forming products and the norm ratio.
  const residual = b.map((v,i) => v/scale - a[i*2]*(x[0]/scale) - a[i*2+1]*(x[1]/scale));
  const denominator = anorm*(xnorm/scale)+bnorm/scale;
  return {backward_error: denominator ? infinity(residual)/denominator : 0,
    residual_infinity: infinity(residual)*scale};
}

export function computeAccuracy(Matrix, config) {
  const a=entries(config.a,4,'A'), b=entries(config.b,2,'b');
  const delta=config.delta, target=config.target;
  if (!Number.isFinite(delta) || Math.abs(delta)>1 || !Number.isInteger(target) || target<0 || target>5)
    throw new RangeError('Choose a valid coefficient and perturbation in [−1, 1]');
  if (!['lu','cholesky','qr','svd'].includes(config.algorithm)) throw new RangeError('Unknown solver');
  if(config.algorithm==='svd' && (!Number.isFinite(config.cutoff)||config.cutoff<0||config.cutoff>1)) throw new RangeError('Relative cutoff must be in [0, 1]');
  const perturbedA=a.slice(), perturbedB=b.slice(), selected=target<4?perturbedA:perturbedB, index=target<4?target:target-4;
  const before=selected[index]; selected[index]+=delta;
  const actualDelta=selected[index]-before;
  function solve(values, rhs) {
    let input, right, factor, output;
    try {
      input=new Matrix(2,2,values); right=new Matrix(2,1,rhs);
      if(config.algorithm==='svd') {
        output=input.solveMinimumNorm(right,config.cutoff);
        const d=input.spectralDiagnostics(config.cutoff), r=input.svd();
        let singular_values;
        try {singular_values=Array.from(r.values);} finally {r.u.dispose();r.vt.dispose();}
        const x=Array.from(output.toArray());
        return {x,rank:d.rank,reciprocal_condition:d.reciprocalCondition,retained_reciprocal_condition:d.retainedReciprocalCondition,
          singular_values,condition_norm:'2-norm',solution_norm:Math.hypot(...x),...diagnostics(values,rhs,x)};
      }
      factor=config.algorithm==='lu'?input.factorLU():config.algorithm==='cholesky'?input.factorCholesky():input.factorQR();
      output=factor.solve(right);
      const x=Array.from(output.toArray());
      return {x, condition_norm:'infinity norm', solution_norm:Math.hypot(...x), reciprocal_condition:factor.reciprocalCondition(), ...diagnostics(values,rhs,x)};
    } finally { output?.dispose(); factor?.dispose(); right?.dispose(); input?.dispose(); }
  }
  const base=solve(a,b);
  let perturbed;
  try {perturbed=solve(perturbedA,perturbedB);} catch(error) {return {a,b,perturbed_a:perturbedA,perturbed_b:perturbedB,base,actual_delta:actualDelta,perturbed_error:error.message};}
  const norm=infinity(base.x), difference=infinity(perturbed.x.map((v,i)=>v-base.x[i]));
  return {a,b,perturbed_a:perturbedA,perturbed_b:perturbedB,base,perturbed,actual_delta:actualDelta,
    relative_solution_change:norm?difference/norm:(difference?null:0)};
}

if (typeof WorkerGlobalScope !== 'undefined' && globalThis instanceof WorkerGlobalScope) {
  let initialization;
  const handle = async ({data}) => {
    const id=data?.id;
    try {
      if(data?.type!=='accuracy') throw new Error('Unsupported accuracy request');
      initialization ||= (async()=>{
        const bundle=data.bundle, urls=[];
        const url=source=>{const value=URL.createObjectURL(new Blob([source],{type:'text/javascript'}));urls.push(value);return value;};
        try {
          const {createMatrixAPI}=await import(url(bundle.wrapper_source));
          const api=await createMatrixAPI({moduleUrl:url(bundle.module_source),wasmBinary:Uint8Array.from(atob(bundle.wasm_base64),c=>c.charCodeAt(0)),locateFile:()=> 'embedded.wasm'});
          const {runChecks}=await import(url(bundle.worker_source));
          globalThis.onmessage=handle;
          const checks=await runChecks(api.Matrix,bundle.fixtures);
          if(checks.passed!==checks.total) throw new Error('Shared correctness checks failed');
          return {...api,checks};
        } finally {globalThis.onmessage=handle; urls.forEach(value=>URL.revokeObjectURL(value));}
      })();
      const {Matrix,checks}=await initialization;
      globalThis.postMessage({type:'accuracy',id,result:computeAccuracy(Matrix,data.config),checks});
    } catch(error) {globalThis.postMessage({type:'error',id,message:error instanceof Error?error.message:String(error)});}
  };
  globalThis.onmessage=handle;
}
