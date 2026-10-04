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
const metricDescriptions={"total": ["Total time", "Shorter bars mean less measured time, including preparation."], "solve": ["Prepared solve time", "Shorter bars mean a faster solve with an existing factor; preparation is excluded."], "setup": ["ILU setup time", "Shorter bars mean faster ILU(0) factor preparation; solving and ordering are excluded."], "iterations": ["GMRES iterations", "Shorter bars mean fewer accepted GMRES iterations; each iteration and preparation still have a time cost."], "bandwidth": ["Structural bandwidth", "Shorter bars mean stored entries stay nearer the diagonal; this does not by itself predict convergence or total time."]};

// Plot coordinates use CSS pixels; resizing never shrinks labels with a fixed viewBox.
function drawSavedChart(orders,values,metric,isCount){
 const svg=$('saved-chart'),width=Math.max(240,svg.clientWidth||700),compact=width<520,row=compact?66:48,height=orders.length*row+8;
 svg.replaceChildren();svg.setAttribute('viewBox',`0 0 ${width} ${height}`);svg.setAttribute('height',height);
 const start=compact?0:78,end=width-(compact?0:132),span=Math.max(0,end-start),max=Math.max(1,...values.filter(Number.isFinite));
 orders.forEach((order,i)=>{
  const top=i*row+4,label=order==='natural'?'Natural':order.toUpperCase(),value=values[i],available=Number.isFinite(value),note=available?(isCount?String(value):ns(value)):(order==='natural'&&metric==='ordering'?'No ordering required':'Unavailable / failed');
  svg.append(svgNode('text',{x:0,y:top+16,fill:'currentColor','font-size':13},label));
  if(available)svg.append(svgNode('rect',{x:start,y:top+(compact?25:0),width:span*value/max,height:22,rx:3,fill:colors[order]}));
  svg.append(svgNode('text',{x:width,y:top+16,fill:'currentColor','font-size':13,'text-anchor':'end'},note));
 });
}
function describeMetric(metric,orders,values,isCount){
 const [label,detail]=metricDescriptions[metric];
 const measured=orders.map((order,i)=>`${order==='natural'?'Natural':order.toUpperCase()}: ${Number.isFinite(values[i])?(isCount?String(values[i]):ns(values[i])):order==='natural'&&metric==='ordering'?'not required':'unavailable'}`).join(' · ');
 $('saved-summary').textContent=`${label}. ${measured}. ${detail}`;
 $('saved-chart').setAttribute('aria-label',`${label}; ${measured}. Exact measurements and deviations are in the table.`);
}
function describeTotal(a,b,label){
 const complete=a?.status==='passed'&&b?.status==='passed'&&Number.isFinite(a.median_ns)&&Number.isFinite(b.median_ns);
 $('total-tradeoff').textContent=complete?`Measured total tradeoff — ${label}: ${ns(b.median_ns)}; Natural: ${ns(a.median_ns)}. ${a.median_ns>0?`${label} / Natural: ${(b.median_ns/a.median_ns).toFixed(2)}×; below 1 means ${label} was faster for this workload.`:'Natural time is zero, so a ratio is not defined.'} This measured end-to-end cost includes preparation; structural improvements alone do not predict it.`:'Measured total tradeoff — a complete total-time comparison is unavailable. Inspect missing or failed measurements in the table.';
}
function canvasFrame(canvas,square=true){
 const width=Math.max(160,canvas.clientWidth||700),height=square?width:Math.max(160,canvas.clientHeight||width*.65),scale=window.devicePixelRatio||1;
 const backingWidth=Math.round(width*scale),backingHeight=Math.round(height*scale);
 if(canvas.width!==backingWidth||canvas.height!==backingHeight){canvas.width=backingWidth;canvas.height=backingHeight;}
 const ctx=canvas.getContext('2d');ctx.setTransform(scale,0,0,scale,0,0);ctx.clearRect(0,0,width,height);
 return {ctx,width,height};
}
function saved(){
 const rows=data.results.filter(r=>r.implementation===$('language').value&&r.problem===$('problem').value&&r.size===Number($('size').value));
 const find=(operation,ordering)=>rows.find(r=>r.operation===operation&&r.ordering===ordering);
 $('saved-values').replaceChildren();
 const tableRow=(label,a,b)=>{const tr=document.createElement('tr');for(const text of [label,a,b]){const td=document.createElement('td');td.textContent=text;tr.append(td);}$('saved-values').append(tr);};
 tableRow('Ordering','Not performed',time(find('rcm','rcm')));tableRow('Matrix permutation','Not performed',time(find('permute','rcm')));
 for(const [key,label] of [['setup','ILU setup'],['solve','Prepared solve'],['total','Measured total']])tableRow(label,time(find(key,'natural')),time(find(key,'rcm')));
 for(const [key,label] of [['solver_iterations','GMRES iterations'],['bandwidth','Structural bandwidth']])tableRow(label,...['natural','rcm'].map(o=>String(find('solve',o)?.[key]??'Unavailable')));
 const metric=$('metric').value,isCount=['iterations','bandwidth'].includes(metric),field=metric==='iterations'?'solver_iterations':'bandwidth';
 const selected=['natural','rcm'].map(o=>find(isCount?'solve':metric,o)),values=selected.map(r=>{const value=r?.status==='passed'?(isCount?r[field]:r.median_ns):null;return Number.isFinite(value)&&value>=0?value:null;});
 drawSavedChart(['natural','rcm'],values,metric,isCount);describeMetric(metric,['natural','rcm'],values,isCount);
 const a=find('total','natural'),b=find('total','rcm');
 describeTotal(a,b,'RCM');
}
for(const id of ['language','problem','size','metric'])$(id).addEventListener('change',saved);saved();
let result,worker,ready=false,request=0,debounce,watchdog,animation,origin=0;
const liveConfig=()=>({size:Number($('live-size').value),contrast:Number($('contrast').value),speed:Number($('speed').value),scrambled:$('live-problem').value==='scrambled',preconditioner:$('live-preconditioner').value});
const ink=()=>getComputedStyle(document.body).color;
function draw(){
 if(!result)return;const t=Number($('progress').value)/1000,c=$('pattern'),{ctx,width:w}=canvasFrame(c),pad=40,span=w-2*pad,n=result.n,cell=span/Math.max(1,n-1);ctx.clearRect(0,0,w,w);ctx.strokeStyle='#879b9180';ctx.strokeRect(pad,pad,span,span);
 for(const [i,j,value] of result.points){const row=i*(1-t)+result.inverse[i]*t,col=j*(1-t)+result.inverse[j]*t;ctx.fillStyle=value<0?colors.rcm:colors.natural;const radius=Math.max(1,Math.min(3,cell*.35));ctx.fillRect(pad+col*cell-radius,pad+row*cell-radius,radius*2,radius*2);}
 ctx.fillStyle=ink();ctx.font='13px system-ui';ctx.fillText(`0`,pad,25);ctx.fillText(`${n-1}`,w-pad-25,25);ctx.fillText(`column →`,w/2-40,25);ctx.save();ctx.translate(20,w/2);ctx.rotate(-Math.PI/2);ctx.fillText('row →',-30,0);ctx.restore();
 const c2=$('convergence'),{ctx:g,width:residualWidth}=canvasFrame(c2),residualPad=48,residualSpan=residualWidth-2*residualPad;const all=[...result.natural.residuals,...result.rcm.residuals,result.natural.threshold,result.rcm.threshold].filter(v=>v>0&&Number.isFinite(v));
 if(!all.length)return;const high=Math.log10(Math.max(...all)),low=Math.min(high-1,Math.log10(Math.min(...all))),steps=Math.max(1,result.natural.residuals.length-1,result.rcm.residuals.length-1);
 const y=v=>residualWidth-residualPad-(Math.log10(Math.max(v,10**low))-low)/(high-low)*residualSpan;
 g.font='12px system-ui';g.fillStyle=ink();g.strokeStyle='#879b9160';for(let i=0;i<=4;i++){const exponent=low+(high-low)*i/4,yy=y(10**exponent);g.beginPath();g.moveTo(residualPad,yy);g.lineTo(residualWidth-residualPad,yy);g.stroke();g.fillText((10**exponent).toExponential(0),2,yy-5);}g.fillText('0',residualPad,residualWidth-12);g.fillText(String(steps),residualWidth-residualPad-25,residualWidth-12);g.fillText('accepted iteration',residualWidth/2-70,residualWidth-12);
 for(const order of ['natural','rcm']){g.strokeStyle=colors[order];g.lineWidth=3;g.beginPath();result[order].residuals.forEach((v,i)=>{const x=residualPad+i/steps*residualSpan,yy=y(v);if(i===0)g.moveTo(x,yy);else g.lineTo(x,yy);});g.stroke();}
 const threshold=result.natural.threshold||result.rcm.threshold;if(threshold>0){g.setLineDash([8,7]);g.strokeStyle=ink();g.lineWidth=1;g.beginPath();g.moveTo(residualPad,y(threshold));g.lineTo(residualWidth-residualPad,y(threshold));g.stroke();g.setLineDash([]);}
}
function stop(){cancelAnimationFrame(animation);animation=undefined;$('play').textContent='Animate ordering';}
function tick(now){const t=(now-origin)/10000;if(t>=1){$('progress').value='1000';draw();if($('loop').checked){origin=now;animation=requestAnimationFrame(tick);}else stop();return;}$('progress').value=String(Math.floor(t*1000));draw();animation=requestAnimationFrame(tick);}
$('play').addEventListener('click',()=>{if(animation!==undefined){stop();return;}if(Number($('progress').value)>=1000)$('progress').value='0';origin=performance.now()-Number($('progress').value)*10;$('play').textContent='Pause';animation=requestAnimationFrame(tick);});
$('progress').addEventListener('input',()=>{stop();draw();});
new MutationObserver(()=>draw()).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme','style']});
function show(r){result=r;$('progress').value='0';$('progress').disabled=false;$('play').disabled=false;$('bandwidth-note').textContent=`${r.n} unknowns · ${r.points.length} stored entries · bandwidth ${r.width} → ${r.newWidth}`;
 const table=document.createElement('table'),body=document.createElement('tbody');table.append(body);const ms=x=>x===undefined?'—':x===0?'Below timer resolution':`${x.toPrecision(3)} ms`;
 for(const [label,a,b] of [['', 'Natural','RCM'],['Outcome',r.natural.reason,r.rcm.reason],['Iterations',r.natural.iterations??'—',r.rcm.iterations??'—'],['Final residual',r.natural.residual?.toExponential(3)??'—',r.rcm.residual?.toExponential(3)??'—'],['ILU setup',ms(r.natural.setupMilliseconds),ms(r.rcm.setupMilliseconds)],['Solve',ms(r.natural.solveMilliseconds),ms(r.rcm.solveMilliseconds)],['Live total',ms(r.natural.totalMilliseconds),ms(r.rcm.totalMilliseconds)]]){const tr=document.createElement('tr');for(const text of [label,a,b]){const td=document.createElement('td');td.textContent=text;tr.append(td);}body.append(tr);}$('live-values').replaceChildren(table);draw();}
function schedule(){clearTimeout(debounce);stop();request++;$('play').disabled=true;$('progress').disabled=true;$('live-state').textContent='Waiting for changes to settle…';const id=request;debounce=setTimeout(()=>{if(!ready)return;$('live-state').textContent='Computing both orderings…';worker.postMessage({type:'solve',id,config:liveConfig()});clearTimeout(watchdog);watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='Live solve timed out. Reload to retry.';},30000);},300);}
for(const id of ['live-size','live-problem','contrast','speed','live-preconditioner'])$(id).addEventListener('change',schedule);
function workerMain(){
 let api,compare,grid;const urls=[];const moduleURL=s=>{const u=URL.createObjectURL(new Blob([s],{type:'text/javascript'}));urls.push(u);return u;};
 self.onmessage=async({data})=>{try{if(data.type==='init'){const b=data.bundle,{createMatrixAPI}=await import(moduleURL(b.wrapper_source));api=await createMatrixAPI({moduleUrl:moduleURL(b.module_source),locateFile:()=> 'embedded.wasm',wasmBinary:Uint8Array.from(atob(b.wasm_base64),c=>c.charCodeAt(0))});const {checkSparse}=await import(moduleURL(b.checks_source));const count=checkSparse(api,b.fixtures);grid=(await import(moduleURL(b.heat_source))).diffusionGrid;compare=(await import(moduleURL(b.ordering_source))).compareOrdering;urls.forEach(u=>URL.revokeObjectURL(u));self.postMessage({type:'ready',count});return;}
 const c=data.config;const result=compare(api,grid(c.size,c.contrast,{diffusivity:.2,speed:c.speed,angle:30}),c);self.postMessage({type:'result',id:data.id,result});}
 catch(e){self.postMessage({type:'error',id:data.id,message:e.message||String(e)});}};
}
if(bundle.available&&bundle.capabilities?.includes('rcm')){
 const url=URL.createObjectURL(new Blob([`(${workerMain.toString()})();`],{type:'text/javascript'}));worker=new Worker(url,{type:'module'});URL.revokeObjectURL(url);
 worker.onmessage=({data:m})=>{if(m.type==='ready'){ready=true;schedule();return;}if(m.id!==undefined&&m.id!==request)return;clearTimeout(watchdog);if(m.type==='error'){$('live-state').textContent=`Live run failed: ${m.message}`;return;}$('live-state').textContent='Both solves finished. Results below include the original-system residual check.';show(m.result);};worker.onerror=e=>{clearTimeout(watchdog);ready=false;$('live-state').textContent=`Live worker failed: ${e.message}`;worker.terminate();};worker.postMessage({type:'init',bundle});watchdog=setTimeout(()=>{worker.terminate();ready=false;$('live-state').textContent='WebAssembly initialization timed out. Reload to retry.';},30000);
}else $('live-state').textContent=bundle.reason||'This report has no verified ordering WebAssembly bundle.';
window.addEventListener('resize',()=>{saved();draw();});
window.addEventListener('pagehide',()=>{worker?.terminate();clearTimeout(debounce);clearTimeout(watchdog);stop();});
