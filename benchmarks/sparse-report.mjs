const $=id=>document.getElementById(id),data=JSON.parse($('data').textContent),bundle=JSON.parse($('live').textContent);
const hasIlu=bundle.capabilities?.includes('ilu0')===true;
const hasGmres=bundle.capabilities?.includes('gmres')||bundle.fixtures?.some(test=>test.op==='gmres');
if(!hasGmres){$('solver').value='cg';$('solver').disabled=true;$('speed').value='0';for(const id of ['speed','angle','restart'])$(id).disabled=true;}
$('ilu-option').disabled=!hasIlu;$('hot').disabled=!hasIlu;$('capability-note').hidden=hasIlu;
$('preset-transport').disabled=!hasGmres;$('preset-ilu').disabled=!hasGmres||!hasIlu;
$('live-preset').value=hasGmres?'transport':'diffusion';
const add=(parent,tag,text)=>{const e=document.createElement(tag);e.textContent=text;parent.append(e);return e;};
const colors=['#31766a','#b95122','#487ab6','#8865a3','#b34f79','#87711f','#2d8292','#8d6046','#6673aa','#986635','#6c7b3d'];
const darkColors=['#82cbb1','#f1ad75','#94b9ed','#c6a2e7','#eda0bb','#d3bc69','#7cc6d1','#d8ad8f','#b1b9ec','#d9af73','#b6c683'];
const names=new Map(data.implementations.map(i=>[i.id,i.name||i.id]));
const markerSymbols=['●','■','◆','▲'];
const dashes=[[],[7,4],[2,4]];
const seriesColor=index=>(document.documentElement.dataset.theme==='dark'?darkColors:colors)[index%colors.length];
const seriesDash=index=>dashes[Math.floor(index/markerSymbols.length)%dashes.length];
let current,frame=0,playing=false,lastTime=0,worker,ready=false,sequence=0,pending,timer,watchdog,workerURL;
let visible=new Set(data.implementations.map(i=>i.id));
const number=x=>Number(x).toLocaleString(undefined,{maximumFractionDigits:6});
const theme=()=>getComputedStyle(document.documentElement);
function canvas(id){const c=$(id),r=c.getBoundingClientRect(),dpr=devicePixelRatio||1;c.width=Math.round(r.width*dpr);c.height=Math.round(r.height*dpr);const ctx=c.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,r.width,r.height);ctx.font='12px system-ui';ctx.fillStyle=theme().getPropertyValue('--muted');return {ctx,w:r.width,h:r.height};}
function drawHeat(){const {ctx,w,h}=canvas('heat');if(!current)return;
 const n=current.size,v=current.iterates[frame],s=Math.min(w-40,h-40),left=(w-s)/2,top=18;
 const stops=[[22,55,104],[52,135,168],[234,214,135],[236,96,43]];
 for(let y=0;y<n;y++)for(let x=0;x<n;x++){const t=Math.max(0,Math.min(1,v[y*n+x]/100))*3,k=Math.min(2,Math.floor(t)),f=t-k;ctx.fillStyle=`rgb(${stops[k].map((a,i)=>Math.round(a+(stops[k+1][i]-a)*f)).join(',')})`;ctx.fillRect(left+x*s/n,top+y*s/n,s/n+.3,s/n+.3);}
 if(current.speed>0){const angle=current.angle*Math.PI/180,cx=w/2,cy=top+s*.7,len=Math.min(30,s*.13),dx=Math.cos(angle),dy=Math.sin(angle);ctx.strokeStyle='rgba(255,255,255,.85)';ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(cx-dx*len,cy-dy*len);ctx.lineTo(cx+dx*len,cy+dy*len);ctx.lineTo(cx+dx*(len-7)+dy*5,cy+dy*(len-7)-dx*5);ctx.moveTo(cx+dx*len,cy+dy*len);ctx.lineTo(cx+dx*(len-7)-dy*5,cy+dy*(len-7)+dx*5);ctx.stroke();ctx.lineWidth=1;}
 ctx.fillStyle=theme().getPropertyValue('--muted');ctx.textAlign='center';ctx.fillText(`${current.hot??100} · hot boundary`,w/2,13);ctx.fillText(`Iteration ${frame} · interior ${n} × ${n}`,w/2,h-8);
}
function drawResidual(){const {ctx,w,h}=canvas('residual');if(!current)return;
 const history=current.residuals,estimates=current.estimatedResiduals||[],positive=[...history,...estimates].filter(x=>x>0),floor=Math.max(Number.MIN_VALUE,Math.min(current.threshold,...positive)*.2),low=Math.log10(floor),high=Math.log10(Math.max(...positive,1)),span=Math.max(1,high-low),left=65,right=w-20,top=20,bottom=h-40;
 const px=i=>left+(right-left)*i/Math.max(1,history.length-1),py=v=>bottom-(Math.log10(Math.max(v,floor))-low)/span*(bottom-top);
 ctx.strokeStyle=theme().getPropertyValue('--line');ctx.fillStyle=theme().getPropertyValue('--muted');ctx.textAlign='right';
 for(let i=0;i<=4;i++){const y=top+(bottom-top)*i/4,value=10**(low+span*(1-i/4));ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText(value.toExponential(0),left-8,y+4);}
 ctx.textAlign='center';ctx.fillText('0',left,bottom+18);ctx.fillText(String(history.length-1),right,bottom+18);ctx.fillText('Solver iteration',w/2,h-6);
 ctx.strokeStyle='#ba8c30';ctx.setLineDash([5,4]);ctx.beginPath();ctx.moveTo(left,py(current.threshold));ctx.lineTo(right,py(current.threshold));ctx.stroke();ctx.setLineDash([]);
 for(const step of current.restarts||[]){ctx.strokeStyle='#ba638b';ctx.setLineDash([2,4]);ctx.beginPath();ctx.moveTo(px(step),top);ctx.lineTo(px(step),bottom);ctx.stroke();}ctx.setLineDash([]);
 if(estimates.length){ctx.strokeStyle='#32a7a0';ctx.setLineDash([6,4]);ctx.beginPath();estimates.forEach((v,i)=>i?ctx.lineTo(px(i),py(v)):ctx.moveTo(px(i),py(v)));ctx.stroke();ctx.setLineDash([]);}
 ctx.strokeStyle=theme().getPropertyValue('--accent');ctx.lineWidth=2;ctx.beginPath();history.forEach((v,i)=>i?ctx.lineTo(px(i),py(v)):ctx.moveTo(px(i),py(v)));ctx.stroke();ctx.beginPath();ctx.arc(px(frame),py(history[frame]),5,0,Math.PI*2);ctx.fillStyle=ctx.strokeStyle;ctx.fill();ctx.lineWidth=1;
}
let iterationViews=[],selectedIterationRow;
const residualText=value=>Number.isFinite(value)?value===0?'0':value.toExponential(3):'Not recorded';
function isGmres(){return current?.solver==='gmres'||!!current?.estimatedResiduals?.length;}
function restartPosition(index){
 if(!isGmres())return 'Not used by CG';
 if(index===0)return 'Initial guess';
 const previous=(current.restarts||[]).filter(step=>step>0&&step<=index).at(-1)||0;
 return previous===index?`Restart at iteration ${index}`:previous?`${index-previous} step${index-previous===1?'':'s'} after restart ${previous}`:`Initial cycle · step ${index}`;
}
function renderIterationTable(){
 $('iteration-rows').replaceChildren();iterationViews=[];selectedIterationRow=undefined;
 for(const [index,value] of current.residuals.entries()){
  const row=add($('iteration-rows'),'tr',''),identity=add(row,'td',''),button=add(identity,'button',String(index));
  button.type='button';button.setAttribute('aria-label',`Show iteration ${index}`);button.disabled=index>=current.iterates.length;
  button.onclick=()=>{setPlaying(false);updateFrame(index);};
  for(const text of [residualText(value),isGmres()?residualText(current.estimatedResiduals?.[index]):'Not applicable',residualText(current.threshold),restartPosition(index)])add(row,'td',text);
  iterationViews.push(row);
 }
}
function updateFrame(value){
 if(!current)return;
 frame=Math.max(0,Math.min(current.iterates.length-1,value));$('frame').value=frame;$('step').textContent=`Iteration ${frame} / ${current.iterations}`;
 $('frame-residual').textContent=residualText(current.residuals[frame]);
 $('frame-estimate').textContent=isGmres()?residualText(current.estimatedResiduals?.[frame]):'Not applicable · CG';
 $('frame-tolerance').textContent=residualText(current.threshold);$('frame-restart').textContent=restartPosition(frame);
 selectedIterationRow?.removeAttribute('aria-current');selectedIterationRow=iterationViews[frame];selectedIterationRow?.setAttribute('aria-current','step');
 drawHeat();drawResidual();
}
function setPlaying(value){playing=value;$('play').textContent=playing?'Pause':'Play';lastTime=0;}
function animate(time){if(playing&&current){if(!lastTime)lastTime=time;const duration=15000/Math.max(1,current.iterations);if(time-lastTime>=duration){const next=frame+Math.max(1,Math.floor((time-lastTime)/duration));lastTime=time;if(next>current.iterations){if($('loop').checked)updateFrame(0);else{updateFrame(current.iterations);setPlaying(false);}}else updateFrame(next);}}requestAnimationFrame(animate);}
requestAnimationFrame(animate);
$('play').onclick=()=>{if(!playing&&current&&frame===current.iterations)updateFrame(0);setPlaying(!playing);};
$('frame').oninput=()=>{setPlaying(false);updateFrame(Number($('frame').value));};
function config(){return {size:Number($('grid').value),contrast:Number($('contrast').value),limit:Number($('limit').value),jacobi:$('preconditioner').value==='jacobi',ilu:$('preconditioner').value==='ilu0',hot:Number($('hot').value),solver:$('solver').value,restart:Number($('restart').value),diffusivity:Number($('diffusivity').value),speed:Number($('speed').value),angle:Number($('angle').value)};}
function fail(message){clearTimeout(watchdog);$('status').textContent=message;setPlaying(false);$('play').disabled=true;$('frame').disabled=true;worker?.terminate();worker=undefined;ready=false;if(workerURL){URL.revokeObjectURL(workerURL);workerURL=undefined;}}
function guard(){clearTimeout(watchdog);watchdog=setTimeout(()=>fail('The solver timed out. Change a setting to restart.'),30000);}
function dispatch(){if(!ready||!pending)return;guard();worker.postMessage(pending);pending=undefined;}
function start(){if(worker)return;if(!bundle.available){fail(bundle.reason);return;}workerURL=URL.createObjectURL(new Blob([bundle.worker_source],{type:'text/javascript'}));worker=new Worker(workerURL,{type:'module'});guard();worker.onerror=e=>fail('Solver failed: '+e.message);worker.onmessage=({data:m})=>{
 if(m.type==='ready'){clearTimeout(watchdog);ready=true;$('verification').textContent=`${m.verified} shared sparse fixtures passed in this browser. Live bundle SHA-256: ${bundle.sha256}`;dispatch();return;}
 if(m.type==='error'&&m.id===undefined){fail('Solver initialization failed: '+m.message);return;}
 if(m.id!==sequence)return;
 clearTimeout(watchdog);if(m.type==='error'){fail('Solver failed: '+m.message);return;}
 current=m.result;$('status').textContent=current.converged?'Converged: the true residual meets the requested tolerance.':`Stopped: ${current.reason.replaceAll('_',' ')}. The displayed estimate has not converged.`;
 $('stats').replaceChildren();for(const [value,label] of [[current.size**2,'unknowns'],[current.nnz,'nonzero coefficients'],[current.iterations,'accepted iterations'],[(current.restarts||[]).length,'restarts'],...(current.ilu?[[current.reusedFactor?'Reused':current.setupMilliseconds>0?number(current.setupMilliseconds)+' ms':'Below timer resolution','ILU setup'],[current.solveMilliseconds>0?number(current.solveMilliseconds)+' ms':'Below timer resolution','live solve']]:[]),[current.residuals.at(-1).toExponential(2),'final true residual']]){const e=add($('stats'),'div','');e.className='stat';add(e,'strong',value);add(e,'span',label);}
 renderIterationTable();$('frame').max=current.iterations;$('frame').disabled=false;$('play').disabled=current.iterations===0;setPlaying(false);updateFrame(current.iterations);
 };worker.postMessage({type:'init',bundle});}
function schedule(){sequence++;pending={type:'solve',id:sequence,config:config()};clearTimeout(timer);setPlaying(false);$('play').disabled=true;$('frame').disabled=true;$('status').textContent='Updating solver…';timer=setTimeout(()=>{start();dispatch();},300);}
const guidedPresets={
 diffusion:{solver:'cg',preconditioner:'jacobi',grid:'16',restart:'20',diffusivity:'0.2',speed:'0',angle:'30',contrast:'1',limit:'400',hot:'100'},
 transport:{solver:'gmres',preconditioner:'jacobi',grid:'16',restart:'20',diffusivity:'0.2',speed:'4',angle:'30',contrast:'1',limit:'400',hot:'100'},
 ilu:{solver:'gmres',preconditioner:'ilu0',grid:'16',restart:'20',diffusivity:'0.2',speed:'4',angle:'30',contrast:'1',limit:'400',hot:'100'}
};
const presetNotes={
 diffusion:'Follow heat diffusion with CG and Jacobi. Scrub the accepted steps to watch the residual change.',
 transport:'Add directed flow and follow restarted GMRES. Compare preconditioners while keeping the problem fixed.',
 ilu:'Build an ILU(0) factor, then change the hot boundary under Problem to reuse it for a new right-hand side.',
 custom:'Custom setup. Adjust the problem or advanced solver settings below the plots.'
};
function syncLiveControls(){
 if(!hasGmres)$('solver').value='cg';
 const cg=$('solver').value==='cg';$('restart').disabled=cg;$('speed').disabled=cg;$('angle').disabled=cg;$('ilu-option').disabled=cg||!hasIlu;
 if(cg)$('speed').value='0';
 if((cg||!hasIlu)&&$('preconditioner').value==='ilu0')$('preconditioner').value='jacobi';
 $('hot-value').textContent=$('hot').value;$('angle-value').textContent=$('angle').value+'°';
 $('preset-note').textContent=presetNotes[$('live-preset').value]||presetNotes.custom;
}
function changed(){ $('live-preset').value='custom';syncLiveControls();schedule(); }
for(const id of ['grid','contrast','limit','preconditioner','hot','solver','restart','diffusivity','speed','angle'])$(id).onchange=changed;
$('angle').oninput=changed;$('hot').oninput=changed;
$('live-preset').onchange=()=>{
 const key=$('live-preset').value;
 if((key==='transport'&&!hasGmres)||(key==='ilu'&&(!hasGmres||!hasIlu))){$('live-preset').value='custom';syncLiveControls();return;}
 const preset=guidedPresets[key];if(!preset){syncLiveControls();return;}
 for(const [id,value] of Object.entries(preset))if(!$(id).disabled||['solver','speed','angle','restart'].includes(id))$(id).value=value;
 syncLiveControls();schedule();
};
syncLiveControls();
function compact(value){
 if(!Number.isFinite(value))return '—';
 if(value===0)return '0';
 return Number(value.toPrecision(3)).toLocaleString(undefined,{maximumFractionDigits:6});
}
function duration(ns){return ns>=1e9?`${compact(ns/1e9)} s`:ns>=1e6?`${compact(ns/1e6)} ms`:ns>=1e3?`${compact(ns/1e3)} µs`:`${compact(ns)} ns`;}
function bytes(value){return value>=1048576?`${compact(value/1048576)} MiB`:value>=1024?`${compact(value/1024)} KiB`:`${number(value)} B`;}
function niceStep(value){
 const magnitude=10**Math.floor(Math.log10(value)),fraction=value/magnitude;
 return (fraction<=1?1:fraction<=2?2:fraction<=2.5?2.5:fraction<=5?5:10)*magnitude;
}
function timingAxis(values,logarithmic){
 const positive=values.filter(v=>v>0&&Number.isFinite(v));
 const maximum=Math.max(...positive,1e-9);
 const unit=maximum>=1e9?{divisor:1e9,label:'s'}:maximum>=1e6?{divisor:1e6,label:'ms'}:maximum>=1e3?{divisor:1e3,label:'µs'}:{divisor:1,label:'ns'};
 if(logarithmic&&positive.length){
  let low=Math.floor(Math.log10(Math.min(...positive))),high=Math.ceil(Math.log10(maximum));
  if(low===high){low--;high++;}
  const step=Math.max(1,Math.ceil((high-low)/5)),ticks=[];
  for(let exponent=low;exponent<=high;exponent+=step)ticks.push(10**exponent);
  if(Math.log10(ticks.at(-1))<high)ticks.push(10**high);
  return {...unit,ticks,position:value=>(Math.log10(value)-low)/(high-low)};
 }
 const step=niceStep(maximum/4),top=Math.ceil(maximum/step)*step,ticks=[];
 for(let index=0;index<=Math.round(top/step);index++)ticks.push(index*step);
 return {...unit,ticks,position:value=>value/top};
}
function axisNumber(value){
 return value!==0&&(Math.abs(value)<.001||Math.abs(value)>=1e6)?value.toExponential(0):compact(value);
}
function marker(ctx,index,x,y,radius){
 ctx.beginPath();
 switch(index%4){
  case 0:ctx.arc(x,y,radius,0,Math.PI*2);break;
  case 1:ctx.rect(x-radius,y-radius,radius*2,radius*2);break;
  case 2:ctx.moveTo(x,y-radius*1.3);ctx.lineTo(x+radius*1.3,y);ctx.lineTo(x,y+radius*1.3);ctx.lineTo(x-radius*1.3,y);ctx.closePath();break;
  default:ctx.moveTo(x,y-radius*1.4);ctx.lineTo(x+radius*1.2,y+radius);ctx.lineTo(x-radius*1.2,y+radius);ctx.closePath();
 }
 ctx.fill();
}
const legendKeys=[],languageInputs=[];
function drawTimings(){
 const op=$('operation').value,rows=data.results.filter(r=>r.operation===op&&visible.has(r.implementation));
 const highlighted=$('highlight').value,logarithmic=$('timing-scale').value==='log';
 $('timing-rows').replaceChildren();
 for(const row of [...rows].sort((a,b)=>(names.get(a.implementation)||a.implementation).localeCompare(names.get(b.implementation)||b.implementation)||a.size-b.size)){
  const tr=add($('timing-rows'),'tr','');
  if(row.implementation===highlighted)tr.className='series-selected';
  for(const value of [names.get(row.implementation)||row.implementation,`${row.size} × ${row.size} / ${number(row.unknowns)}`,row.status==='passed'?duration(row.median_ns):row.status,row.status==='passed'?duration(row.mad_ns):'—',row.solver_iterations??row.cg_iterations??'—',row.logical_workspace_bytes===undefined?'—':bytes(row.logical_workspace_bytes),row.logical_preconditioner_bytes===undefined?'—':bytes(row.logical_preconditioner_bytes)])add(tr,'td',value);
 }
 const {ctx,w,h}=canvas('timings'),valid=rows.filter(r=>r.status==='passed'&&Number.isFinite(r.median_ns));
 const plotted=valid.filter(r=>!logarithmic||r.median_ns>0);
 const maxX=Math.max(1,...plotted.map(r=>r.unknowns)),axis=timingAxis(plotted.map(r=>r.median_ns),logarithmic);
 const left=66,right=w-24,top=30,bottom=h-48;
 const px=value=>left+(right-left)*value/maxX,py=value=>bottom-(bottom-top)*axis.position(value);
 ctx.strokeStyle=theme().getPropertyValue('--line');ctx.fillStyle=theme().getPropertyValue('--muted');ctx.textAlign='right';
 for(const tick of axis.ticks){const y=py(tick);ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText(axisNumber(tick/axis.divisor),left-9,y+4);}
 ctx.fillText(axis.label,left-9,16);ctx.textAlign='center';
 const sizes=[...new Set(plotted.map(r=>r.unknowns))].sort((a,b)=>a-b);
 const stride=Math.max(1,Math.ceil(sizes.length/Math.max(2,Math.floor((right-left)/75))));
 sizes.forEach((n,index)=>{if(index%stride===0||index===sizes.length-1)ctx.fillText(number(n),px(n),bottom+20);});
 ctx.fillText('Number of unknowns · grid width²',w/2,h-8);
 // Draw the selected series last so crossings do not conceal it.
 const series=[...data.implementations.entries()].sort((a,b)=>(a[1].id===highlighted)-(b[1].id===highlighted));
 for(const [index,impl] of series){
  const points=plotted.filter(r=>r.implementation===impl.id).sort((a,b)=>a.unknowns-b.unknowns);
  ctx.globalAlpha=highlighted&&highlighted!==impl.id ? .2 : 1;
  ctx.strokeStyle=ctx.fillStyle=seriesColor(index);ctx.lineWidth=impl.id===highlighted?3:2;ctx.setLineDash(seriesDash(index));ctx.beginPath();
  points.forEach((row,i)=>i?ctx.lineTo(px(row.unknowns),py(row.median_ns)):ctx.moveTo(px(row.unknowns),py(row.median_ns)));ctx.stroke();ctx.setLineDash([]);
  for(const row of points)marker(ctx,index,px(row.unknowns),py(row.median_ns),impl.id===highlighted?5:3.5);
 }
 ctx.globalAlpha=1;
 for(const [index,key] of legendKeys.entries())key.style.color=seriesColor(index);
 if(!plotted.length){ctx.fillStyle=theme().getPropertyValue('--muted');ctx.fillText('No included measurements',w/2,h/2);}
 const omitted=valid.length-plotted.length;
 $('timing-note').textContent=`${visible.size}/${data.implementations.length} languages shown · ${logarithmic?'Logarithmic time axis':'Linear time axis from zero'} · scale fits visible measurements${highlighted?' · Highlighting '+names.get(highlighted):''}${omitted?' · '+omitted+' zero timing(s) shown only in the table':''}.`;
}
for(const [index,impl] of data.implementations.entries()){
 const option=add($('highlight'),'option',impl.name||impl.id);option.value=impl.id;
 const label=add($('languages'),'label',''),input=document.createElement('input');input.type='checkbox';input.checked=true;languageInputs.push(input);
 input.onchange=()=>{input.checked?visible.add(impl.id):visible.delete(impl.id);if(!visible.has($('highlight').value))$('highlight').value='';drawTimings();};
 const key=document.createElement('span');key.className='series-key '+['','dashed','dotted'][Math.floor(index/4)%3];key.setAttribute('aria-hidden','true');add(key,'span',markerSymbols[index%4]);legendKeys.push(key);
 label.append(input,key,document.createTextNode(impl.name||impl.id));
}
const workloads=[...new Set(data.results.map(r=>r.operation))];
if(workloads.length){$('operation').replaceChildren();for(const op of workloads){const option=add($('operation'),'option',({ilu_setup:'ILU(0) · setup',ilu_apply:'ILU(0) · triangular solves',gmres_none:'GMRES · no preconditioner',gmres_jacobi:'GMRES · Jacobi',gmres_ilu_reused:'GMRES · saved ILU(0)',gmres_ilu_total:'GMRES · ILU(0) setup + solve'}[op]||op.replaceAll('_',' ').replace('gmres','GMRES').replace('jacobi','+ Jacobi').replace('csr spmv','CSR matrix-vector product').replace('dense spmv','Dense matrix-vector product').replace('cg','CG')));option.value=op;}$('operation').value=workloads.includes('csr_spmv')?'csr_spmv':workloads[0];}
$('operation').onchange=drawTimings;$('timing-scale').onchange=drawTimings;
$('highlight').onchange=()=>{
 const id=$('highlight').value;
 if(id&&!visible.has(id)){
  visible.add(id);const index=data.implementations.findIndex(impl=>impl.id===id);languageInputs[index].checked=true;
 }
 drawTimings();
};
const storage=new Map(data.results.map(r=>[r.size,r]));
for(const row of [...storage.values()].sort((a,b)=>a.size-b.size)){
 const tr=add($('storage'),'tr','');
 for(const value of [`${row.size} × ${row.size} / ${number(row.unknowns)}`,number(row.nnz),bytes(row.logical_dense_bytes),bytes(row.logical_csr_bytes),compact(row.logical_dense_bytes/row.logical_csr_bytes)+'×'])add(tr,'td',value);
}
const measured=data.results.filter(row=>row.status==='passed'),gridCount=new Set(data.results.map(row=>row.size)).size;
for(const [value,label] of [[`${data.implementations.filter(i=>i.status==='passed').length}/${data.implementations.length}`,'implementations verified'],[number(measured.length),'timed workloads'],[gridCount,'grid sizes measured']]){
 const card=add($('overview'),'div','');add(card,'strong',value);add(card,'span',label);
}
$('run-context').textContent=`Recorded ${data.created_at?.slice(0,10)||'date unavailable'} · ${data.machine.os||'OS unspecified'} ${data.machine.architecture||''} · revision ${data.revision?.slice(0,12)||'unavailable'}${data.dirty?' · working changes':''}`;
$('provenance').textContent=JSON.stringify({created_at:data.created_at,machine:data.machine,revision:data.revision,uncommitted_changes:data.dirty,source_sha256:data.source_sha256},null,2);
for(const [key,value] of Object.entries(data.methodology)){
 const item=add($('method'),'p','');add(item,'strong',key.replaceAll('_',' ').replace(/^./,c=>c.toUpperCase()));item.append(document.createTextNode(value));
}
for(const impl of data.implementations){
 add($('toolchains'),'h3',`${impl.name||impl.id} · ${impl.status} · ${impl.checks?.filter(c=>c.passed).length||0} checks`);add($('toolchains'),'p',impl.toolchain||'');
 for(const command of impl.build_commands||[])add($('toolchains'),'pre',command.join(' '));
}

function drawReuse(){
 if(data.suite!=='ilu-reuse-v1')return;
 const id=$('reuse-language').value,size=Number($('reuse-grid').value),count=Number($('rhs-count').value)||8;$('rhs-value').textContent=count;
 const rows=new Map(data.results.filter(r=>r.implementation===id&&r.size===size&&r.status==='passed').map(r=>[r.operation,r.median_ns]));
 const {ctx,w,h}=canvas('reuse-chart'),setup=rows.get('ilu_setup'),solve=rows.get('gmres_ilu_reused');
 if(setup===undefined||solve===undefined){$('reuse-values').replaceChildren();$('reuse-note').textContent='No complete measurements for this selection.';return;}
 const series=[['None',rows.get('gmres_none')*count],['Jacobi',rows.get('gmres_jacobi')*count],['ILU rebuilt',rows.get('gmres_ilu_total')*count],['ILU reused',setup+solve*count]].filter(([,value])=>Number.isFinite(value));
 $('reuse-values').replaceChildren();for(const [label,value] of series)add($('reuse-values'),'li',`${label}: ${duration(value)}`);
 const maximum=Math.max(1,...series.map(s=>s[1])),compactLayout=w<600,left=compactLayout?16:180,right=w-(compactLayout?16:100);
 series.forEach(([label,value],i)=>{
  const y=24+i*(compactLayout?66:55);ctx.fillStyle=theme().getPropertyValue('--muted');
  if(compactLayout){ctx.textAlign='left';ctx.fillText(label,left,y,w*.5);ctx.textAlign='right';ctx.fillText(duration(value),right,y,w*.4);}
  else{ctx.textAlign='right';ctx.fillText(label,left-10,y+15);ctx.textAlign='left';ctx.fillText(duration(value),left+(right-left)*value/maximum+8,y+16);}
  ctx.fillStyle=seriesColor(i);ctx.fillRect(left,y+(compactLayout?10:0),Math.max(1,(right-left)*value/maximum),24);
 });
 const saving=rows.get('gmres_jacobi')-solve,breakEven=saving>0?Math.floor(setup/saving)+1:null;
 $('reuse-note').textContent=`${count} RHS: setup + ${count} × reused solve. ${breakEven?`Estimated to beat Jacobi after ${breakEven} RHS.`:'These medians do not predict an advantage over Jacobi.'} Setup ${duration(setup)}; reused solve ${duration(solve)}. GC, warmup and workload changes can alter the result.`;
}
if(data.suite==='ilu-reuse-v1'){
 $('reuse-section').hidden=false;$('reuse-nav').hidden=false;
 for(const impl of data.implementations){const option=add($('reuse-language'),'option',impl.name||impl.id);option.value=impl.id;}
 $('reuse-language').value=data.implementations[0]?.id||'';
 const widths=[...new Set(data.results.map(r=>r.size))].sort((a,b)=>a-b);for(const width of widths){const option=add($('reuse-grid'),'option',String(width));option.value=String(width);}$('reuse-grid').value=String(widths[0]);
 for(const id of ['reuse-language','reuse-grid','rhs-count'])$(id).onchange=drawReuse;$('rhs-count').oninput=drawReuse;
}

function redraw(){drawHeat();drawResidual();drawTimings();drawReuse();}
window.addEventListener('matrix-theme-change',redraw);
new ResizeObserver(redraw).observe(document.querySelector('main'));
window.addEventListener('pagehide',()=>{clearTimeout(timer);clearTimeout(watchdog);worker?.terminate();if(workerURL)URL.revokeObjectURL(workerURL);});
redraw();schedule();
