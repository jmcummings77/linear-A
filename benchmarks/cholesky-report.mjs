const $=id=>document.getElementById(id),data=JSON.parse($('data').textContent),bundle=JSON.parse($('live').textContent);
const names=Object.fromEntries(data.implementations.map(i=>[i.id,i.name]));
const option=(select,value,label=value)=>{const o=document.createElement('option');o.value=value;o.textContent=label;select.append(o);};
Object.entries(names).forEach(([id,name])=>option($('language'),id,name));
[...new Set(data.results.map(r=>r.size))].sort((a,b)=>a-b).forEach(n=>option($('size'),n));
$('provenance').textContent=`${data.created_at} · ${data.machine.os} ${data.machine.architecture} · revision ${data.revision.slice(0,12)}${data.dirty?' · working changes':''} · ${data.implementations.filter(i=>i.status==='passed').length}/${data.implementations.length} ports verified`;
for(const [key,value] of Object.entries(data.methodology)){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=value;$('method').append(dt,dd);}
$('toolchains').textContent=JSON.stringify(data.implementations,null,2);
const ns=x=>x<1e3?`${x.toPrecision(3)} ns`:x<1e6?`${(x/1e3).toPrecision(3)} µs`:`${(x/1e6).toPrecision(3)} ms`;
const time=row=>!row?'Not measured':row.status!=='passed'?`Failed: ${row.error||row.status}`:`${ns(row.median_ns)} ± ${ns(row.mad_ns)}`;
const colors={natural:'#b79240',rcm:'#479b8b'};
function svgNode(tag,attrs,text){const e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,String(v));if(text!==undefined)e.textContent=text;return e;}
function saved(){
 const rows=data.results.filter(r=>r.implementation===$('language').value&&r.problem===$('problem').value&&r.size===Number($('size').value));
 const find=(operation,ordering)=>rows.find(r=>r.operation===operation&&r.ordering===ordering);
 $('saved-values').replaceChildren();
 const tableRow=(label,a,b)=>{const tr=document.createElement('tr');for(const text of [label,a,b]){const td=document.createElement('td');td.textContent=text;tr.append(td);}$('saved-values').append(tr);};
 for(const [key,label] of [['symbolic','Symbolic analysis'],['factor','Numerical factorization'],['solve','Prepared solve'],['total','Measured total']])tableRow(label,time(find(key,'natural')),time(find(key,'rcm')));
 for(const [key,label] of [['fill_count','Fill entries'],['factor_nnz','Lower-factor entries'],['logical_factor_bytes','Logical factor bytes']])tableRow(label,...['natural','rcm'].map(o=>String(find('solve',o)?.[key]??'Unavailable')));
 const metric=$('metric').value,isCount=['fill_count','factor_nnz','logical_factor_bytes'].includes(metric),field=metric;
 const selected=['natural','rcm'].map(o=>find(isCount?'solve':metric,o)),values=selected.map(r=>r?.status==='passed'?(isCount?r[field]:r.median_ns):null),max=Math.max(1,...values.filter(v=>v!==null));
 const svg=$('saved-chart');svg.replaceChildren();svg.setAttribute('viewBox','0 0 1000 260');
 ['natural','rcm'].forEach((o,i)=>{const y=65+i*95;svg.append(svgNode('text',{x:0,y:y+20,fill:'currentColor'},o==='natural'?'Natural':'RCM'));if(values[i]!==null){const width=620*values[i]/max;svg.append(svgNode('rect',{x:100,y,width,height:35,rx:4,fill:colors[o]}));svg.append(svgNode('text',{x:115+width,y:y+23,fill:'currentColor'},isCount?values[i]:ns(values[i])));}else svg.append(svgNode('text',{x:100,y:y+23,fill:'currentColor'},'Unavailable / failed'));});
 const a=find('total','natural'),b=find('total','rcm');
 $('saved-summary').textContent=a?.status==='passed'&&b?.status==='passed'?`RCM total / natural total: ${(b.median_ns/a.median_ns).toFixed(2)}×. Below 1 means RCM was faster in this recorded workload. Less fill alone does not predict this ratio.`:'A complete total-time comparison is not available; inspect the failures above.';
}
for(const id of ['language','problem','size','metric'])$(id).addEventListener('change',saved);saved();
let result,worker,ready=false,request=0,debounce,watchdog,animation,origin=0;
const liveConfig=()=>({size:Number($('live-size').value),contrast:Number($('contrast').value),scrambled:$('live-problem').value==='scrambled'});
const ink=()=>getComputedStyle(document.body).color;
function draw(){
 if(!result)return;const step=Math.min(result.n,Math.floor(Number($('progress').value)/1000*result.n));
 const counts=[];
 for(const [order,id] of [['natural','pattern'],['rcm','convergence']]){
  const r=result[order],c=$(id),ctx=c.getContext('2d'),w=c.width,pad=40,span=w-2*pad,cell=span/Math.max(1,result.n-1);ctx.clearRect(0,0,w,w);ctx.strokeStyle='#879b9180';ctx.strokeRect(pad,pad,span,span);
  let count=0,fresh=0;
  for(const [i,j,birth] of r.points){if(birth>=step)continue;if(birth>=0)count++;if(birth===step-1&&birth>=0)fresh++;
   ctx.fillStyle=birth===step-1&&birth>=0?'#e17b62':birth<0?colors.natural:colors.rcm;
   const radius=Math.max(1,Math.min(4,cell*.35));ctx.fillRect(pad+j*cell-radius,pad+i*cell-radius,radius*2,radius*2);
  }
  ctx.fillStyle=ink();ctx.font='16px system-ui';ctx.fillText(`0`,pad,25);ctx.fillText(`${result.n-1}`,w-pad-25,25);ctx.fillText('column →',w/2-40,25);
  if(step>0){const at=pad+(step-1)*cell;ctx.strokeStyle='#e17b6280';ctx.beginPath();ctx.moveTo(pad,at);ctx.lineTo(w-pad,at);ctx.moveTo(at,pad);ctx.lineTo(at,w-pad);ctx.stroke();}
  counts.push(`${order==='rcm'?'RCM':'Natural'}: ${count} fill entries (${fresh} new)`);
 }
 $('bandwidth-note').textContent=`${result.n} unknowns · ${step===0?'Before elimination':`After pivot ${step-1}`} · ${counts.join(' · ')}`;
 $('frame-code').textContent=step===0?'Original lower structure + diagonal; no variables eliminated.':`pivot k = ${step-1}\nremaining = neighbors(k) with index > k\nfor each pair (i, j) in remaining: add edge (i, j) if absent\nrecord new entries with fillStep = ${step-1}`;
}

function stop(){cancelAnimationFrame(animation);animation=undefined;$('play').textContent='Animate elimination';}
function tick(now){const t=(now-origin)/10000;if(t>=1){$('progress').value='1000';draw();if($('loop').checked){origin=now;animation=requestAnimationFrame(tick);}else stop();return;}$('progress').value=String(Math.floor(t*1000));draw();animation=requestAnimationFrame(tick);}
$('play').addEventListener('click',()=>{if(animation!==undefined){stop();return;}if(Number($('progress').value)>=1000)$('progress').value='0';origin=performance.now()-Number($('progress').value)*10;$('play').textContent='Pause';animation=requestAnimationFrame(tick);});
$('progress').addEventListener('input',()=>{stop();draw();});
new MutationObserver(()=>draw()).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme','style']});
function show(r){result=r;$('progress').value='0';for(const id of ['progress','play','previous','next'])$(id).disabled=false;
 const table=document.createElement('table'),body=document.createElement('tbody');table.append(body);const ms=x=>x===undefined?'—':x===0?'Below timer resolution':`${x.toPrecision(3)} ms`;
 const rows=[['','Natural','RCM'],['Outcome',...['natural','rcm'].map(o=>r[o].ok?'Verified':r[o].error)]];
 for(const [field,label] of [['nnz','Lower-factor entries'],['fillCount','Fill entries'],['logicalFactorBytes','Logical factor bytes']])rows.push([label,...['natural','rcm'].map(o=>r[o][field]??'—')]);
 for(const [field,label] of [['residual','Relative solution residual'],['reconstructionError','Relative reconstruction error']])rows.push([label,...['natural','rcm'].map(o=>r[o][field]?.toExponential(3)??'—')]);
 for(const [field,label] of [['orderingMilliseconds','Ordering + permutations'],['symbolicMilliseconds','Symbolic analysis'],['factorMilliseconds','Numerical factorization'],['solveMilliseconds','Prepared solve'],['totalMilliseconds','Live total']])rows.push([label,...['natural','rcm'].map(o=>ms(r[o][field]))]);
 for(const row of rows){const tr=document.createElement('tr');for(const text of row){const td=document.createElement('td');td.textContent=text;tr.append(td);}body.append(tr);}$('live-values').replaceChildren(table);draw();
}
for(const [id,delta] of [['previous',-1],['next',1]])$(id).addEventListener('click',()=>{if(!result)return;stop();const step=Math.floor(Number($('progress').value)*result.n/1000);$('progress').value=String(Math.ceil(Math.max(0,Math.min(result.n,step+delta))*1000/result.n));draw();});

function schedule(){clearTimeout(debounce);stop();request++;for(const id of ['play','progress','previous','next'])$(id).disabled=true;$('live-state').textContent='Waiting for changes to settle…';const id=request;debounce=setTimeout(()=>{if(!ready)return;$('live-state').textContent='Computing both orderings…';worker.postMessage({type:'solve',id,config:liveConfig()});clearTimeout(watchdog);watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='Live solve timed out. Reload to retry.';},30000);},300);}
for(const id of ['live-size','live-problem','contrast'])$(id).addEventListener('change',schedule);
function workerMain(){
 let api,compare,grid;const urls=[];const moduleURL=s=>{const u=URL.createObjectURL(new Blob([s],{type:'text/javascript'}));urls.push(u);return u;};
 self.onmessage=async({data})=>{try{if(data.type==='init'){const b=data.bundle,{createMatrixAPI}=await import(moduleURL(b.wrapper_source));api=await createMatrixAPI({moduleUrl:moduleURL(b.module_source),locateFile:()=> 'embedded.wasm',wasmBinary:Uint8Array.from(atob(b.wasm_base64),c=>c.charCodeAt(0))});const {checkSparse}=await import(moduleURL(b.checks_source));const count=checkSparse(api,b.fixtures);grid=(await import(moduleURL(b.heat_source))).diffusionGrid;compare=(await import(moduleURL(b.cholesky_source))).compareCholesky;urls.forEach(u=>URL.revokeObjectURL(u));self.postMessage({type:'ready',count});return;}
 const c=data.config;const result=compare(api,grid(c.size,c.contrast,{diffusivity:1,speed:0}),c);self.postMessage({type:'result',id:data.id,result});}
 catch(e){self.postMessage({type:'error',id:data.id,message:e.message||String(e)});}};
}
if(bundle.available&&bundle.capabilities?.includes('cholesky')){
 const url=URL.createObjectURL(new Blob([`(${workerMain.toString()})();`],{type:'text/javascript'}));worker=new Worker(url,{type:'module'});URL.revokeObjectURL(url);
 worker.onmessage=({data:m})=>{if(m.type==='ready'){ready=true;schedule();return;}if(m.id!==undefined&&m.id!==request)return;clearTimeout(watchdog);if(m.type==='error'){$('live-state').textContent=`Live run failed: ${m.message}`;return;}$('live-state').textContent='Both factorizations finished. Inspect verification outcomes below.';show(m.result);};worker.onerror=e=>{clearTimeout(watchdog);ready=false;$('live-state').textContent=`Live worker failed: ${e.message}`;worker.terminate();};worker.postMessage({type:'init',bundle});watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='WebAssembly initialization timed out. Reload to retry.';},30000);
}else $('live-state').textContent=bundle.reason||'This report has no verified Cholesky WebAssembly bundle.';
window.addEventListener('pagehide',()=>{worker?.terminate();clearTimeout(debounce);clearTimeout(watchdog);stop();});
