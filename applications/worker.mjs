let ready;
globalThis.onmessage=async({data})=>{
  const {id,kind,config,bundle}=data;
  try{
    ready ||= (async()=>{
      const urls=[],url=s=>{const u=URL.createObjectURL(new Blob([s],{type:'text/javascript'}));urls.push(u);return u;};
      try{
        const {createMatrixAPI}=await import(url(bundle.wrapper_source));
        const api=await createMatrixAPI({moduleUrl:url(bundle.module_source),wasmBinary:Uint8Array.from(atob(bundle.wasm_base64),c=>c.charCodeAt(0)),locateFile:()=> 'embedded.wasm'});
        const math=await import(url(bundle.math_source));return {...api,math};
      }finally{urls.forEach(u=>URL.revokeObjectURL(u));}
    })();
    const {Matrix,math}=await ready;
    const result=kind==='pca'?math.pca(Matrix,config.points):kind==='fit'?math.fit(Matrix,config.points,config.degree):kind==='image'?math.compress(Matrix,config.pixels,config.size,config.rank):(()=>{throw new Error('Unknown application');})();
    postMessage({id,kind,result});
  }catch(error){postMessage({id,kind,error:error.message});}
};
