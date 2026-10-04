const $=id=>document.getElementById(id),data=JSON.parse($('data').textContent),bundle=JSON.parse($('live').textContent);
const names=Object.fromEntries(data.implementations.filter(i=>i.id!=='suitesparse').map(i=>[i.id,i.name]));
const option=(select,value,label=value)=>{const o=document.createElement('option');o.value=value;o.textContent=label;select.append(o);};
Object.entries(names).forEach(([id,name])=>option($('language'),id,name));
[...new Set(data.results.map(r=>r.size))].sort((a,b)=>a-b).forEach(n=>option($('size'),n));
$('provenance').textContent=`${data.created_at} · ${data.machine.os} ${data.machine.architecture} · revision ${data.revision.slice(0,12)}${data.dirty?' · working changes':''} · ${data.implementations.filter(i=>i.status==='passed'&&i.id!=='suitesparse').length}/${data.implementations.filter(i=>i.id!=='suitesparse').length} ports verified`;
for(const [key,value] of Object.entries(data.methodology)){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=value;$('method').append(dt,dd);}
$('toolchains').textContent=JSON.stringify(data.implementations,null,2);
const ns=x=>x<1e3?`${x.toPrecision(3)} ns`:x<1e6?`${(x/1e3).toPrecision(3)} µs`:`${(x/1e6).toPrecision(3)} ms`;
const time=row=>!row?'Not measured':row.status!=='passed'?`Failed: ${row.error||row.status}`:`${ns(row.median_ns)} ± ${ns(row.mad_ns)}`;
const colors={natural:'#b79240',rcm:'#479b8b',amd:'#9273bf'};
function svgNode(tag,attrs,text){const e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,String(v));if(text!==undefined)e.textContent=text;return e;}
function saved(){
 const rows=data.results.filter(r=>r.implementation===$('language').value&&r.problem===$('problem').value&&r.size===Number($('size').value));
 const find=(operation,ordering)=>rows.find(r=>r.operation===operation&&r.ordering===ordering);
 $('saved-values').replaceChildren();
 const tableRow=(label,...values)=>{const tr=document.createElement('tr');for(const text of [label,...values]){const td=document.createElement('td');td.textContent=text;tr.append(td);}$('saved-values').append(tr);};
 for(const [key,label] of [['ordering','Ordering only'],['symbolic','Symbolic analysis'],['factor','Numerical factorization'],['solve','Prepared solve'],['total','Measured total']])tableRow(label,...['natural','rcm','amd'].map(o=>o==='natural'&&key==='ordering'?'Not needed':time(find(key,o))));
 for(const [key,label] of [['fill_count','Fill entries'],['factor_nnz','Lower-factor entries'],['logical_factor_bytes','Logical factor bytes']])tableRow(label,...['natural','rcm','amd'].map(o=>String(find('solve',o)?.[key]??'Unavailable')));
 const metric=$('metric').value,isCount=['fill_count','factor_nnz','logical_factor_bytes'].includes(metric),field=metric;
 const selected=['natural','rcm','amd'].map(o=>find(isCount?'solve':metric,o)),values=selected.map(r=>r?.status==='passed'?(isCount?r[field]:r.median_ns):null),max=Math.max(1,...values.filter(v=>v!==null));
 const svg=$('saved-chart');svg.replaceChildren();svg.setAttribute('viewBox','0 0 1000 350');
 ['natural','rcm','amd'].forEach((o,i)=>{const y=65+i*95;svg.append(svgNode('text',{x:0,y:y+20,fill:'currentColor'},o==='natural'?'Natural':o.toUpperCase()));if(values[i]!==null){const width=620*values[i]/max;svg.append(svgNode('rect',{x:100,y,width,height:35,rx:4,fill:colors[o]}));svg.append(svgNode('text',{x:115+width,y:y+23,fill:'currentColor'},isCount?values[i]:ns(values[i])));}else svg.append(svgNode('text',{x:100,y:y+23,fill:'currentColor'},o==='natural'&&metric==='ordering'?'No ordering required':'Unavailable / failed'));});
 const a=find('total','natural'),b=find('total','amd');
 $('saved-summary').textContent=a?.status==='passed'&&b?.status==='passed'?`AMD total / natural total: ${(b.median_ns/a.median_ns).toFixed(2)}×. Below 1 means AMD was faster in this recorded workload. Less fill alone does not predict this ratio.`:'A complete total-time comparison is not available; inspect the failures above.';
 const ref=data.results.find(r=>r.implementation==='suitesparse'&&r.problem===$('problem').value&&r.size===Number($('size').value));
 $('reference').textContent=ref?`SuiteSparse AMD reference: ${ref.fill_count} fill entries · ${ref.factor_nnz} lower-factor entries · ${ref.logical_factor_bytes} logical bytes · ordering ${time(ref)}. Different heuristics and ties; exact permutation equality is not expected. Reference timing includes Python ctypes overhead; its samples run separately in one process.`:'SuiteSparse reference was not requested for this report.';
}
for(const id of ['language','problem','size','metric'])$(id).addEventListener('change',saved);saved();
let result,worker,ready=false,request=0,debounce,watchdog,animation,origin=0;
const liveConfig=()=>({size:Number($('live-size').value),contrast:Number($('contrast').value),scrambled:$('live-problem').value==='scrambled',problem:$('live-problem').value});
const ink=()=>getComputedStyle(document.body).color;
function draw(){
 if(!result)return;const step=Math.min(result.n,Math.floor(Number($('progress').value)/1000*result.n));
 const counts=[];
 for(const [order,id] of [['natural','pattern'],['rcm','convergence'],['amd','amd-pattern']]){
  const r=result[order],c=$(id),ctx=c.getContext('2d'),w=c.width,pad=40,span=w-2*pad,cell=span/Math.max(1,result.n-1);ctx.clearRect(0,0,w,w);ctx.strokeStyle='#879b9180';ctx.strokeRect(pad,pad,span,span);
  let count=0,fresh=0;
  for(const [i,j,birth] of r.points){if(birth>=step)continue;if(birth>=0)count++;if(birth===step-1&&birth>=0)fresh++;
   ctx.fillStyle=birth===step-1&&birth>=0?'#e17b62':birth<0?colors.natural:colors.rcm;
   const radius=Math.max(1,Math.min(4,cell*.35));ctx.fillRect(pad+j*cell-radius,pad+i*cell-radius,radius*2,radius*2);
  }
  ctx.fillStyle=ink();ctx.font='16px system-ui';ctx.fillText(`0`,pad,25);ctx.fillText(`${result.n-1}`,w-pad-25,25);ctx.fillText('column →',w/2-40,25);
  if(step>0){const at=pad+(step-1)*cell;ctx.strokeStyle='#e17b6280';ctx.beginPath();ctx.moveTo(pad,at);ctx.lineTo(w-pad,at);ctx.moveTo(at,pad);ctx.lineTo(at,w-pad);ctx.stroke();}
  counts.push(`${order==='natural'?'Natural':order.toUpperCase()}: ${count} fill entries (${fresh} new)`);
 }
 $('bandwidth-note').textContent=`${result.n} unknowns · ${step===0?'Before elimination':`After pivot ${step-1}`} · ${counts.join(' · ')}`;
 drawGraph(step);

}

function stop(){cancelAnimationFrame(animation);animation=undefined;$('play').textContent='Animate elimination';}
function tick(now){const t=(now-origin)/10000;if(t>=1){$('progress').value='1000';draw();if($('loop').checked){origin=now;animation=requestAnimationFrame(tick);}else stop();return;}$('progress').value=String(Math.floor(t*1000));draw();animation=requestAnimationFrame(tick);}
$('play').addEventListener('click',()=>{if(animation!==undefined){stop();return;}if(Number($('progress').value)>=1000)$('progress').value='0';origin=performance.now()-Number($('progress').value)*10;$('play').textContent='Pause';animation=requestAnimationFrame(tick);});
$('progress').addEventListener('input',()=>{stop();draw();});
new MutationObserver(()=>draw()).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme','style']});
function show(r){result=r;$('progress').value='0';for(const id of ['progress','play','previous','next'])$(id).disabled=false;
 const table=document.createElement('table'),body=document.createElement('tbody');table.append(body);const ms=x=>x===undefined?'—':x===0?'Below timer resolution':`${x.toPrecision(3)} ms`;
 const rows=[['','Natural','RCM','AMD'],['Outcome',...['natural','rcm','amd'].map(o=>r[o].ok?'Verified':r[o].error)]];
 for(const [field,label] of [['nnz','Lower-factor entries'],['fillCount','Fill entries'],['logicalFactorBytes','Logical factor bytes']])rows.push([label,...['natural','rcm','amd'].map(o=>r[o][field]??'—')]);
 for(const [field,label] of [['residual','Relative solution residual'],['reconstructionError','Relative reconstruction error']])rows.push([label,...['natural','rcm','amd'].map(o=>r[o][field]?.toExponential(3)??'—')]);
 for(const [field,label] of [['orderingMilliseconds','Ordering + permutations'],['symbolicMilliseconds','Symbolic analysis'],['factorMilliseconds','Numerical factorization'],['solveMilliseconds','Prepared solve'],['totalMilliseconds','Live total']])rows.push([label,...['natural','rcm','amd'].map(o=>ms(r[o][field]))]);
 for(const row of rows){const tr=document.createElement('tr');for(const text of row){const td=document.createElement('td');td.textContent=text;tr.append(td);}body.append(tr);}$('live-values').replaceChildren(table);draw();
}
for(const [id,delta] of [['previous',-1],['next',1]])$(id).addEventListener('click',()=>{if(!result)return;stop();const step=Math.floor(Number($('progress').value)*result.n/1000);$('progress').value=String(Math.ceil(Math.max(0,Math.min(result.n,step+delta))*1000/result.n));draw();});

function schedule(){clearTimeout(debounce);stop();request++;for(const id of ['play','progress','previous','next'])$(id).disabled=true;$('live-state').textContent='Waiting for changes to settle…';const id=request;debounce=setTimeout(()=>{if(!ready)return;$('live-state').textContent='Computing three orderings…';worker.postMessage({type:'solve',id,config:liveConfig()});clearTimeout(watchdog);watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='Live solve timed out. Reload to retry.';},30000);},300);}
for(const id of ['live-size','live-problem','contrast'])$(id).addEventListener('change',schedule);
function workerMain(){
 let api,compare,grid,graph;const urls=[];const moduleURL=s=>{const u=URL.createObjectURL(new Blob([s],{type:'text/javascript'}));urls.push(u);return u;};
 self.onmessage=async({data})=>{try{if(data.type==='init'){const b=data.bundle,{createMatrixAPI}=await import(moduleURL(b.wrapper_source));api=await createMatrixAPI({moduleUrl:moduleURL(b.module_source),locateFile:()=> 'embedded.wasm',wasmBinary:Uint8Array.from(atob(b.wasm_base64),c=>c.charCodeAt(0))});const {checkSparse}=await import(moduleURL(b.checks_source));const count=checkSparse(api,b.fixtures);grid=(await import(moduleURL(b.heat_source))).diffusionGrid;const live=await import(moduleURL(b.amd_source));compare=live.compareAMD;graph=live.graphProblem;urls.forEach(u=>URL.revokeObjectURL(u));self.postMessage({type:'ready',count});return;}
 const c=data.config;const result=compare(api,['tree','irregular'].includes(c.problem)?graph(c.size,c.problem):grid(c.size,c.contrast,{diffusivity:1,speed:0}),c);self.postMessage({type:'result',id:data.id,result});}
 catch(e){self.postMessage({type:'error',id:data.id,message:e.message||String(e)});}};
}
if(bundle.available&&bundle.capabilities?.includes('amd')){
 const url=URL.createObjectURL(new Blob([`(${workerMain.toString()})();`],{type:'text/javascript'}));worker=new Worker(url,{type:'module'});URL.revokeObjectURL(url);
 worker.onmessage=({data:m})=>{if(m.type==='ready'){ready=true;schedule();return;}if(m.id!==undefined&&m.id!==request)return;clearTimeout(watchdog);if(m.type==='error'){$('live-state').textContent=`Live run failed: ${m.message}`;return;}$('live-state').textContent='All three factorizations finished. Inspect verification outcomes below.';show(m.result);};worker.onerror=e=>{clearTimeout(watchdog);ready=false;$('live-state').textContent=`Live worker failed: ${e.message}`;worker.terminate();};worker.postMessage({type:'init',bundle});watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='WebAssembly initialization timed out. Reload to retry.';},30000);
}else $('live-state').textContent=bundle.reason||'This report has no verified AMD WebAssembly bundle.';
window.addEventListener('pagehide',()=>{worker?.terminate();clearTimeout(debounce);clearTimeout(watchdog);stop();});

function drawGraph(step){
 const kind=$('graph-order').value,r=result[kind],canvas=$('graph'),ctx=canvas.getContext('2d'),w=canvas.width,h=canvas.height,n=result.n;
 ctx.clearRect(0,0,w,h);
 if(!r.ok){$('frame-code').textContent=r.error;return;}
 const positions=Array.from({length:n},(_,i)=>{const t=2*Math.PI*i/n-Math.PI/2;return [w/2+(w/2-50)*Math.cos(t),h/2+(h/2-50)*Math.sin(t)];});
 const eliminated=new Set(r.order.slice(0,Math.max(0,step-1))),frame=step?r.trace.frames[step-1]:null,pivot=frame?.pivot,neighbors=new Set(frame?.neighbors||[]);
 for(const [i,j,birth] of r.trace.edges){
  if(birth>=step)continue;const past=eliminated.has(i)||eliminated.has(j),fresh=birth>=0&&birth===step-1;
  if(past&&!$('past-edges').checked)continue;
  ctx.strokeStyle=past?'#879b9120':fresh?'#e17b62':birth>=0?'#479b8b99':'#b7924060';ctx.lineWidth=fresh?3:1;
  ctx.beginPath();ctx.moveTo(...positions[i]);ctx.lineTo(...positions[j]);ctx.stroke();
 }
 for(let i=0;i<n;i++){
  ctx.fillStyle=i===pivot?'#e17b62':eliminated.has(i)?'#879b9140':neighbors.has(i)?'#479b8b':colors[kind];ctx.beginPath();ctx.arc(...positions[i],i===pivot?10:neighbors.has(i)?6:4,0,Math.PI*2);ctx.fill();
  if(n<=64){ctx.fillStyle=ink();ctx.font='12px system-ui';ctx.fillText(String(i),positions[i][0]+7,positions[i][1]-7);}
 }
 $('graph-note').textContent=`${kind==='natural'?'Natural':kind.toUpperCase()} · Original vertex labels and fixed positions. Coral pivot / new edges; teal remaining neighbors / earlier fill; faint eliminated vertices.`;
 $('frame-code').textContent=frame?`Step ${step} / ${n}: eliminate original vertex ${pivot}
Remaining neighbors: [${frame.neighbors.join(', ')}]
New fill edges (${frame.fill.length}): ${frame.fill.map(e=>`(${e.join(', ')})`).join(' ')||'none'}
Connect remaining neighbors, then remove the pivot.
WASM factor pattern agrees with independent explicit-graph replay.`:'No vertices eliminated. Choose an ordering, then step or animate to follow actual fill edges.';
}
for(const id of ['graph-order','past-edges'])$(id).addEventListener('change',draw);
