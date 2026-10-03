/** Steady-state diffusion, with solver iterates (not physical time steps). */
export function diffusionGrid(size, contrast=0) {
  if (!Number.isInteger(size) || size < 2 || size > 32 || !Number.isFinite(contrast) || contrast < 0 || contrast > 3) throw new RangeError('invalid diffusion grid');
  const n=size*size, k=Array.from({length:n},(_,i)=>10**(contrast*(i%size)/(size-1)));
  const offsets=[0], indices=[], values=[], b=Array(n).fill(0);
  for(let i=0;i<n;i++) {
    const y=Math.floor(i/size),x=i%size, entries=[];let diagonal=0;
    for(const [dy,dx] of [[-1,0],[0,-1],[0,1],[1,0]]) {
      const yy=y+dy,xx=x+dx;let weight;
      if(yy>=0&&yy<size&&xx>=0&&xx<size) {const j=yy*size+xx;weight=Math.sqrt(k[i]*k[j]);entries.push([j,-weight]);}
      else {weight=k[i];if(yy<0)b[i]+=100*weight;}
      diagonal+=weight;
    }
    entries.push([i,diagonal]);entries.sort((a,b)=>a[0]-b[0]);
    for(const [j,v] of entries){indices.push(j);values.push(v);}offsets.push(values.length);
  }
  return {size,n,offsets,indices,values,b};
}
export function solveHeat(api,{size=16,contrast=2,jacobi=true,limit=400}={}) {
  if(typeof jacobi!=='boolean'||!Number.isInteger(limit)||limit<0||limit>800)throw new RangeError('invalid solver settings');
  const grid=diffusionGrid(size,contrast);
  const a=new api.CSRMatrix(grid.n,grid.n,grid.offsets,grid.indices,grid.values);
  try {
    const result=a.conjugateGradient(grid.b,{rtol:1e-8,atol:0,maxIterations:limit,jacobi,capture:true});
    return {...result,size,nnz:grid.values.length,threshold:1e-8*Math.hypot(...grid.b)};
  } finally {a.dispose?.();}
}
