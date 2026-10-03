/** Steady advection–diffusion; frames are solver iterates, not physical time. */
export function diffusionGrid(size, contrast=0, {diffusivity=1,speed=0,angle=0}={}) {
  if (!Number.isInteger(size) || size < 2 || size > 32 || !Number.isFinite(contrast) || contrast < 0 || contrast > 3 || !Number.isFinite(diffusivity)||diffusivity<.01||diffusivity>2||!Number.isFinite(speed)||speed<0||speed>10||!Number.isFinite(angle)||angle< -180||angle>180) throw new RangeError('invalid transport grid');
  const n=size*size, k=Array.from({length:n},(_,i)=>10**(contrast*(i%size)/(size-1)));
  const offsets=[0], indices=[], values=[], b=Array(n).fill(0),theta=angle*Math.PI/180,vx=speed*Math.cos(theta)/(size+1),vy=speed*Math.sin(theta)/(size+1);
  for(let i=0;i<n;i++) {
    const y=Math.floor(i/size),x=i%size, entries=new Map();let diagonal=0;
    for(const [dy,dx] of [[-1,0],[0,-1],[0,1],[1,0]]) {
      const yy=y+dy,xx=x+dx;let weight;
      if(yy>=0&&yy<size&&xx>=0&&xx<size) {const j=yy*size+xx;weight=diffusivity*Math.sqrt(k[i]*k[j]);entries.set(j,-weight);}
      else {weight=diffusivity*k[i];if(yy<0)b[i]+=100*weight;}
      diagonal+=weight;
    }
    // First-order upwinding for v·grad(T), after multiplying the PDE by h².
    for(const [component,dy,dx] of [[vx,0,vx>=0?-1:1],[vy,vy>=0?-1:1,0]]) {
      const weight=Math.abs(component),yy=y+dy,xx=x+dx;diagonal+=weight;
      if(yy>=0&&yy<size&&xx>=0&&xx<size){const j=yy*size+xx;entries.set(j,(entries.get(j)||0)-weight);}
      else if(yy<0)b[i]+=100*weight;
    }
    entries.set(i,diagonal);
    for(const [j,v] of [...entries].sort((a,b)=>a[0]-b[0])){indices.push(j);values.push(v);}offsets.push(values.length);
  }
  return {size,n,offsets,indices,values,b};
}
export function solveHeat(api,{size=16,contrast=2,jacobi=true,limit=400,solver='cg',restart=20,diffusivity=1,speed=0,angle=30,ilu=false,hot=100}={},cache) {
  if(typeof jacobi!=='boolean'||!Number.isInteger(limit)||limit<0||limit>800||!['cg','gmres'].includes(solver)||!Number.isInteger(restart)||restart<1||restart>100)throw new RangeError('invalid solver settings');
  if(typeof ilu!=='boolean'||!Number.isFinite(hot)||hot<0||hot>100||(ilu&&(jacobi||solver!=='gmres')))throw new RangeError('invalid preconditioner or boundary');
  if(solver==='cg'&&speed!==0)throw new RangeError('CG requires symmetric diffusion: set flow speed to zero or select GMRES.');
  const grid=diffusionGrid(size,contrast,{diffusivity,speed,angle});
  const a=new api.CSRMatrix(grid.n,grid.n,grid.offsets,grid.indices,grid.values);
  let factor,setupMilliseconds=0,reusedFactor=false;
  try {
    if(ilu){
      const key=JSON.stringify([size,contrast,diffusivity,speed,angle]);
      if(cache?.key===key&&cache.factor){factor=cache.factor;reusedFactor=true;}
      else{const start=performance.now();factor=new api.ILU0(a);setupMilliseconds=performance.now()-start;if(cache){cache.factor?.dispose?.();cache.factor=factor;cache.key=key;}}
    }
    const b=grid.b.map(v=>v*hot/100);
    const options={restart,rtol:1e-8,atol:0,maxIterations:limit,jacobi,capture:true,preconditioner:factor};
    const start=performance.now();
    const result=solver==='gmres'?a.gmres(b,options):a.conjugateGradient(b,options);
    return {...result,setupMilliseconds,solveMilliseconds:performance.now()-start,reusedFactor,ilu,hot,estimatedResiduals:result.estimatedResiduals||[],restarts:result.restarts||[],solver,restart,speed,angle,diffusivity,size,nnz:grid.values.length,threshold:1e-8*Math.hypot(...b)};
  } finally {if(!cache)factor?.dispose?.();a.dispose?.();}
}

/** Cache owned factors while only RHS, restart, tolerance or playback settings change. */
export function createHeatSolver(api){const cache={};const solve=config=>solveHeat(api,config,cache);solve.dispose=()=>{cache.factor?.dispose?.();delete cache.factor;delete cache.key;};return solve;}
