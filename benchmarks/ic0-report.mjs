const $=id=>document.getElementById(id),data=JSON.parse($('data').textContent),bundle=JSON.parse($('live').textContent);
const labels={ic0_setup:'IC(0) setup + export',ic0_apply:'Prepared IC(0) apply',cg:'Plain CG',jacobi:'Jacobi-CG',ic0_solve:'Prepared IC(0)-CG',ic0_total:'IC(0)-CG total',cholesky_solve:'Prepared Cholesky',cholesky_total:'Cholesky total'};
const methods={cg:'Plain CG',jacobi:'Jacobi-CG',ic0:'IC(0)-CG',cholesky:'Full Cholesky'};
const option=(select,value,label=value)=>{const o=document.createElement('option');o.value=value;o.textContent=label;select.append(o);};
data.implementations.forEach(i=>option($('language'),i.id,i.name));[...new Set(data.results.map(r=>r.size))].forEach(n=>option($('size'),n));
$('provenance').textContent=`${data.created_at} · ${data.machine.os} ${data.machine.architecture} · ${data.revision.slice(0,12)}${data.dirty?' · working changes':''} · ${data.implementations.filter(i=>i.status==='passed').length}/${data.implementations.length} ports verified`;
for(const [key,value] of Object.entries(data.methodology)){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=value;$('method').append(dt,dd);}$('toolchains').textContent=JSON.stringify(data.implementations,null,2);
const ns=x=>x<1e3?`${x.toPrecision(3)} ns`:x<1e6?`${(x/1e3).toPrecision(3)} µs`:`${(x/1e6).toPrecision(3)} ms`;
function row(parent,values,header=false){const tr=document.createElement('tr');for(const [index,value] of values.entries()){const cell=document.createElement(header||index===0?'th':'td');if(header||index===0)cell.setAttribute('scope',header?'col':'row');cell.textContent=value;tr.append(cell);}parent.append(tr);}
function svgNode(tag,attrs,text){const el=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))el.setAttribute(k,v);if(text!==undefined)el.textContent=text;return el;}
function saved(){
 const rows=data.results.filter(r=>r.implementation===$('language').value&&r.size===Number($('size').value)&&r.ordering===$('ordering').value);
 const successful=rows.filter(r=>r.status==='passed'),max=Math.max(1,...successful.map(r=>r.median_ns));
 $('saved-values').replaceChildren();
 const chart=$('saved-chart'),width=Math.max(160,chart.getBoundingClientRect().width||600),compact=width<640;
 const entries=Object.entries(labels),rowHeight=compact?66:46,height=entries.length*rowHeight+20;
 const left=compact?0:195,right=width-(compact?0:94),span=Math.max(1,right-left);
 chart.replaceChildren();chart.setAttribute('viewBox',`0 0 ${width} ${height}`);chart.style.height=height+'px';
 for(const [index,[operation,label]] of entries.entries()){
  const measurement=rows.find(r=>r.operation===operation),ok=measurement?.status==='passed',y=12+index*rowHeight;
  row($('saved-values'),[label,ok?`${ns(measurement.median_ns)} ± ${ns(measurement.mad_ns)}`:measurement?.error||'Unavailable',measurement?.solver_iterations??'—',measurement?.factor_nnz??'—']);
  const text=(x,at,value,anchor='start')=>svgNode('text',{x,y:at,fill:'currentColor','font-size':13,'font-family':'system-ui, sans-serif','text-anchor':anchor},value);
  chart.append(text(0,y+(compact?12:18),label));
  if(ok){
   const barWidth=span*measurement.median_ns/max,barY=y+(compact?23:0);
   chart.append(svgNode('rect',{x:left,y:barY,width:span,height:22,rx:3,fill:'var(--track)'}),svgNode('rect',{x:left,y:barY,width:barWidth,height:22,rx:3,fill:'var(--green)'}));
   chart.append(text(compact?width:right+10,y+(compact?12:18),ns(measurement.median_ns),compact?'end':'start'));
  }else chart.append(text(compact?0:left,y+(compact?39:18),'Unavailable / failed'));
 }
 $('saved-summary').textContent=`${successful.length}/${entries.length} operations measured. Median per operation; lower is faster. Bars start at zero and share a scale for this selection. Failures and variation are in the table.`;
}
for(const id of ['language','size','ordering'])$(id).addEventListener('change',saved);saved();
let result,worker,pending=true,ready=false,request=0,debounce,watchdog,animation,origin=0,solutionScale=0,residualScale=0;
const ink=()=>getComputedStyle(document.body).color;
function stop(){cancelAnimationFrame(animation);animation=undefined;$('play').textContent='Play';}
function field(id,values,scale){const c=$(id),ctx=c.getContext('2d'),cell=c.width/result.size;ctx.clearRect(0,0,c.width,c.height);values.forEach((v,i)=>{const t=Math.max(-1,Math.min(1,v/Math.max(scale,1e-300))),a=Math.abs(t),target=t<0?[58,122,170]:[220,85,65];ctx.fillStyle=`rgb(${target.map((c,k)=>Math.round([245,240,219][k]*(1-a)+c*a)).join(',')})`;ctx.fillRect((i%result.size)*cell,Math.floor(i/result.size)*cell,cell+.2,cell+.2);});}
// The preserved IC(0) worker starts at zero and uses rtol=1e-9, atol=0.
const relativeTolerance=1e-9;
function tolerance(){return Math.max(0,...Object.values(result.runs).filter(r=>r.ok).map(r=>r.residuals[0]))*relativeTolerance;}
function curveColor(method){
 const dark=document.documentElement.dataset.theme==='dark';
 return (dark?{cg:'#dfc371',jacobi:'#8fc8de',ic0:'#83cdb9'}:{cg:'#98711f',jacobi:'#36788d',ic0:'#297b69'})[method];
}
const curveDashes={cg:[],jacobi:[7,4],ic0:[2,4]};
function curveMarker(ctx,method,x,y,radius=3){
 ctx.beginPath();
 if(method==='cg')ctx.arc(x,y,radius,0,Math.PI*2);
 else if(method==='jacobi')ctx.rect(x-radius,y-radius,radius*2,radius*2);
 else{ctx.moveTo(x,y-radius*1.3);ctx.lineTo(x+radius*1.3,y);ctx.lineTo(x,y+radius*1.3);ctx.lineTo(x-radius*1.3,y);ctx.closePath();}
 ctx.fill();
}
function curveLayout(){
 const canvas=$('convergence'),width=Math.max(160,canvas.getBoundingClientRect().width||600),height=width<600?300:320,dpr=globalThis.devicePixelRatio||1;
 canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);canvas.style.height=height+'px';
 const ctx=canvas.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,width,height);ctx.font='12px system-ui';
 return {ctx,width,height};
}
function curves(){
 const {ctx,width,height}=curveLayout(),runs=Object.entries(result.runs).filter(([method,r])=>r.ok&&method!=='cholesky');
 const threshold=tolerance(),positive=runs.flatMap(([,r])=>r.residuals).filter(value=>Number.isFinite(value)&&value>0);
 if(threshold>0)positive.push(threshold);
 const hasZero=runs.some(([,r])=>r.residuals.includes(0)),maxStep=Math.max(1,...runs.map(([,r])=>r.iterations));
 const left=64,right=width-18,top=28,bottom=height-(hasZero?76:56),zeroY=bottom+20;
 let low=positive.length?Math.floor(Math.log10(Math.min(...positive))):0,high=positive.length?Math.ceil(Math.log10(Math.max(...positive))):1;
 if(high===low){low--;high++;}
 const px=index=>left+(right-left)*index/maxStep,pyLog=exponent=>bottom-(exponent-low)/(high-low)*(bottom-top),py=value=>pyLog(Math.log10(value));
 const step=Math.max(1,Math.ceil((high-low)/5)),exponents=[];
 for(let exponent=low;exponent<=high;exponent+=step)exponents.push(exponent);
 if(exponents.at(-1)!==high)exponents.push(high);
 ctx.fillStyle=ink();ctx.textAlign='right';ctx.strokeStyle=getComputedStyle(document.documentElement).getPropertyValue('--line');ctx.lineWidth=1;
 for(const exponent of exponents){const y=pyLog(exponent);ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText('1e'+exponent,left-9,y+4);}
 ctx.textAlign='left';ctx.fillText('True residual norm',left,15);
 if(hasZero){ctx.setLineDash([2,4]);ctx.beginPath();ctx.moveTo(left,zeroY);ctx.lineTo(right,zeroY);ctx.stroke();ctx.setLineDash([]);ctx.textAlign='right';ctx.fillText('zero',left-9,zeroY+4);}
 ctx.textAlign='center';const tickCount=width<450?2:4;
 for(let i=0;i<=tickCount;i++){const value=Math.round(maxStep*i/tickCount);if(i>0&&value===Math.round(maxStep*(i-1)/tickCount))continue;ctx.fillText(String(value),px(value),height-26);}
 ctx.fillText('Accepted updates',left+(right-left)/2,height-7);
 if(threshold>0){ctx.strokeStyle=ink();ctx.lineWidth=1.5;ctx.setLineDash([9,5]);ctx.beginPath();ctx.moveTo(left,py(threshold));ctx.lineTo(right,py(threshold));ctx.stroke();ctx.setLineDash([]);}
 const selected=$('method-select').value;
 for(const [method,r] of [...runs].sort((a,b)=>(a[0]===selected)-(b[0]===selected))){
  ctx.strokeStyle=ctx.fillStyle=curveColor(method);ctx.lineWidth=method===selected?2.5:1.5;ctx.setLineDash(curveDashes[method]);ctx.beginPath();let connected=false;
  r.residuals.forEach((value,index)=>{if(value<=0){connected=false;return;}if(connected)ctx.lineTo(px(index),py(value));else ctx.moveTo(px(index),py(value));connected=true;});ctx.stroke();ctx.setLineDash([]);
  const stride=Math.max(1,Math.ceil(r.residuals.length/8));
  r.residuals.forEach((value,index)=>{if(value===0||index%stride===0||index===r.residuals.length-1)curveMarker(ctx,method,px(index),value===0?zeroY:py(value));});
  if(method===selected){const index=Math.min(r.residuals.length-1,Number($('progress').value)),value=r.residuals[index];curveMarker(ctx,method,px(index),value===0?zeroY:py(value),5);}
 }
 const unavailable=Object.entries(result.runs).filter(([method,r])=>method!=='cholesky'&&!r.ok).map(([method])=>methods[method]);
 $('convergence-note').textContent=`Logarithmic axis fits all positive recorded residuals and the stopping tolerance (${threshold.toExponential(3)} = 10⁻⁹ × initial residual). ${hasZero?'Exact zeros use the separate zero row; the log curves stop there. ':''}${selected!=='cholesky'&&result.runs[selected]?.ok?'Larger markers show the selected frame. ':''}Full Cholesky has no CG iteration curve.${unavailable.length?' Unavailable: '+unavailable.join(', ')+'.':''}`;
}
function updateScaleLegend(){
 for(const [id,scale] of [['solution',solutionScale],['residual',residualScale]]){
  $(id+'-min').textContent='−'+scale.toPrecision(3);$(id+'-max').textContent='+'+scale.toPrecision(3);
  $(id+'-scale-note').textContent='Blue: negative · cream: zero · red: positive. Fixed scale across methods and frames.';
 }
}
function draw(){
 if(!result)return;
 const method=$('method-select').value,r=result.runs[method];curves();updateScaleLegend();
 $('current-threshold').textContent=tolerance().toExponential(4);
 if(!r.ok){
  for(const id of ['solution','residual']){const canvas=$(id);canvas.getContext('2d').clearRect(0,0,canvas.width,canvas.height);}
  $('frame-note').textContent=`${methods[method]}: no captured frames`;$('current-residual').textContent='Unavailable';$('completed-status').textContent='Failed: '+r.error;
  $('frame-code').textContent='No factor or solution frames are available for this method.';return;
 }
 const k=Math.min(r.frames.length-1,Number($('progress').value));
 field('solution',r.frames[k],solutionScale);field('residual',r.residualFields[k],residualScale);
 $('frame-note').textContent=`${methods[method]} · ${method==='cholesky'?(k===0?'Initial field':'Final field'):'Accepted update '+k} · frame ${k}/${r.frames.length-1}`;
 $('current-residual').textContent=r.residuals[k].toExponential(4);
 $('completed-status').textContent=`${r.reason}${method==='cholesky'?' · direct solve':' after '+r.iterations+' updates'}`;
 $('frame-code').textContent=method==='cholesky'?'Cholesky: L y = b; Lᵀ x = y. Initial and final fields only.':`Accepted update ${k}\nr = b − A x\nz = ${method==='ic0'?'L⁻ᵀ L⁻¹ r':method==='jacobi'?'diag(A)⁻¹ r':'r'}\nρ = rᵀ z; α = ρ / (pᵀ A p); x ← x + α p\nRecompute the true residual and check convergence before continuing.`;
}
function select(){stop();const r=result?.runs[$('method-select').value],enabled=!pending&&!!r?.ok;$('progress').max=String(enabled?r.frames.length-1:0);$('progress').value='0';for(const id of ['progress','play','previous','next'])$(id).disabled=!enabled;draw();}
function show(r){pending=false;result=r;
 const scaleOf=key=>Object.values(result.runs).filter(run=>run.ok).reduce((scale,run)=>run[key].reduce((value,frame)=>frame.reduce((largest,entry)=>Math.max(largest,Math.abs(entry)),value),scale),0);
 solutionScale=scaleOf('frames');residualScale=scaleOf('residualFields');
 $('live-values').replaceChildren();const table=document.createElement('table');row(table,['Method','Setup (ms)','Solve (ms)','Iterations / status','Factor entries'],true);for(const [method,r] of Object.entries(result.runs))row(table,[methods[method],r.setupMs.toPrecision(3),r.ok?r.solveMs.toPrecision(3):'—',r.ok?`${method==='cholesky'?'direct':r.iterations} / ${r.reason}`:`failed: ${r.error}`,r.nnz??'—']);$('live-values').append(table);$('live-state').textContent=`${r.n} unknowns · ${r.breakdown?'Breakdown example ready.':'Independent true-residual checks passed.'}`;$('live-context').textContent=`Verified ${r.n}-unknown system · ${r.ordering} ordering and coordinate preparation ${r.orderingMs.toPrecision(3)} ms.${r.breakdown?' Positive definite does not guarantee IC(0) success.':''}`;select();}
function tick(now){const t=(now-origin)/10000;$('progress').value=String(Math.min(Number($('progress').max),Math.floor(t*Number($('progress').max))));draw();if(t>=1){if($('loop').checked){origin=now;animation=requestAnimationFrame(tick);}else stop();}else animation=requestAnimationFrame(tick);}
$('play').addEventListener('click',()=>{if(animation!==undefined){stop();return;}if(Number($('progress').value)>=Number($('progress').max))$('progress').value='0';origin=performance.now()-Number($('progress').value)/Math.max(1,Number($('progress').max))*10000;$('play').textContent='Pause';animation=requestAnimationFrame(tick);});$('progress').addEventListener('input',()=>{stop();draw();});$('method-select').addEventListener('change',select);for(const [id,delta] of [['previous',-1],['next',1]])$(id).addEventListener('click',()=>{stop();$('progress').value=String(Math.max(0,Math.min(Number($('progress').max),Number($('progress').value)+delta)));draw();});
function schedule(){pending=true;clearTimeout(debounce);stop();request++;for(const id of ['play','progress','previous','next'])$(id).disabled=true;$('live-state').textContent='Waiting for changes to settle…';const id=request;debounce=setTimeout(()=>{if(!ready)return;$('live-state').textContent='Solving four methods in WebAssembly…';worker.postMessage({type:'solve',id,config:{size:Number($('live-size').value),contrast:Number($('contrast').value),ordering:$('live-ordering').value,breakdown:$('live-problem').value==='breakdown'}});clearTimeout(watchdog);watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='Live solve timed out. Reload to retry.';},30000);},300);}
for(const id of ['live-size','contrast','live-ordering','live-problem'])$(id).addEventListener('change',schedule);
function workerMain(){let api,compare,grid;const urls=[],url=s=>{const u=URL.createObjectURL(new Blob([s],{type:'text/javascript'}));urls.push(u);return u;};self.onmessage=async({data})=>{try{if(data.type==='init'){const b=data.bundle,{createMatrixAPI}=await import(url(b.wrapper_source));api=await createMatrixAPI({moduleUrl:url(b.module_source),locateFile:()=> 'embedded.wasm',wasmBinary:Uint8Array.from(atob(b.wasm_base64),c=>c.charCodeAt(0))});const {checkSparse}=await import(url(b.checks_source));checkSparse(api,b.fixtures);grid=(await import(url(b.heat_source))).diffusionGrid;compare=(await import(url(b.ic0_source))).compareIC0;urls.forEach(u=>URL.revokeObjectURL(u));self.postMessage({type:'ready'});return;}const c=data.config;self.postMessage({type:'result',id:data.id,result:compare(api,grid(c.size,c.contrast),c)});}catch(e){self.postMessage({type:'error',id:data.id,message:e.message||String(e)});}};}
if(bundle.available&&bundle.capabilities?.includes('ic0')){const url=URL.createObjectURL(new Blob([`(${workerMain.toString()})();`],{type:'text/javascript'}));worker=new Worker(url,{type:'module'});URL.revokeObjectURL(url);worker.onmessage=({data:m})=>{if(m.type==='ready'){ready=true;schedule();return;}if(m.id!==undefined&&m.id!==request)return;clearTimeout(watchdog);if(m.type==='error'){$('live-state').textContent=`Live run failed: ${m.message}`;return;}show(m.result);};worker.onerror=e=>{clearTimeout(watchdog);ready=false;$('live-state').textContent=`Live worker failed: ${e.message}`;worker.terminate();};worker.postMessage({type:'init',bundle});watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='WebAssembly initialization timed out. Reload to retry.';},30000);}else $('live-state').textContent=bundle.reason||'No verified IC(0) WebAssembly bundle.';
window.addEventListener('pagehide',()=>{worker?.terminate();clearTimeout(debounce);clearTimeout(watchdog);stop();});

window.addEventListener('matrix-theme-change',()=>draw());
window.addEventListener('resize',()=>{saved();draw();});
