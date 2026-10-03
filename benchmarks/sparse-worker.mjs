let api,solveHeat,verified;
const urls=[];
const moduleURL=source=>{const url=URL.createObjectURL(new Blob([source],{type:'text/javascript'}));urls.push(url);return url;};
self.onmessage=async({data})=>{
  try {
    if(data.type==='init') {
      const bundle=data.bundle;
      const moduleUrl=moduleURL(bundle.module_source);
      const {createMatrixAPI}=await import(moduleURL(bundle.wrapper_source));
      api=await createMatrixAPI({moduleUrl,locateFile:()=> 'embedded.wasm',wasmBinary:Uint8Array.from(atob(bundle.wasm_base64),c=>c.charCodeAt(0))});
      const {checkSparse}=await import(moduleURL(bundle.checks_source));
      verified=checkSparse(api,bundle.fixtures);
      const {createHeatSolver}=await import(moduleURL(bundle.heat_source));solveHeat=createHeatSolver(api);
      for(const url of urls)URL.revokeObjectURL(url);urls.length=0;
      self.postMessage({type:'ready',verified});return;
    }
    if(!api||!solveHeat)throw new Error('solver not initialized');
    self.postMessage({type:'result',id:data.id,result:solveHeat(data.config),verified});
  } catch(error) {self.postMessage({type:'error',id:data.id,message:error.message||String(error)});}
};
