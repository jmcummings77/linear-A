const $=id=>document.getElementById(id),data=JSON.parse($('data').textContent),bundle=JSON.parse($('live').textContent);
const add=(parent,tag,text)=>{const e=document.createElement(tag);e.textContent=text;parent.append(e);return e;};
const colors=['#4682e8','#e87138','#55a878','#b276d2','#dd5d83','#b79724','#439fa8','#936c4d','#7382ad','#be643f','#898d38'];
let current,frame=0,playing=false,lastTime=0,worker,ready=false,sequence=0,pending,timer,watchdog,workerURL;
let visible=new Set(data.implementations.map(i=>i.id));
const number=x=>Number(x).toLocaleString(undefined,{maximumFractionDigits:6});
const theme=()=>getComputedStyle(document.documentElement);
function canvas(id){const c=$(id),r=c.getBoundingClientRect(),dpr=devicePixelRatio||1;c.width=Math.round(r.width*dpr);c.height=Math.round(r.height*dpr);const ctx=c.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,r.width,r.height);ctx.font='12px system-ui';ctx.fillStyle=theme().getPropertyValue('--muted');return {ctx,w:r.width,h:r.height};}
function drawHeat(){const {ctx,w,h}=canvas('heat');if(!current)return;
 const n=current.size,v=current.iterates[frame],s=Math.min(w-40,h-40),left=(w-s)/2,top=18;
 const stops=[[22,55,104],[52,135,168],[234,214,135],[236,96,43]];
 for(let y=0;y<n;y++)for(let x=0;x<n;x++){const t=Math.max(0,Math.min(1,v[y*n+x]/100))*3,k=Math.min(2,Math.floor(t)),f=t-k;ctx.fillStyle=`rgb(${stops[k].map((a,i)=>Math.round(a+(stops[k+1][i]-a)*f)).join(',')})`;ctx.fillRect(left+x*s/n,top+y*s/n,s/n+.3,s/n+.3);}
 ctx.fillStyle=theme().getPropertyValue('--muted');ctx.textAlign='center';ctx.fillText('100 · hot boundary',w/2,13);ctx.fillText(`Iteration ${frame} · interior ${n} × ${n}`,w/2,h-8);
}
function drawResidual(){const {ctx,w,h}=canvas('residual');if(!current)return;
 const history=current.residuals,positive=history.filter(x=>x>0),floor=Math.max(Number.MIN_VALUE,Math.min(current.threshold,...positive)*.2),low=Math.log10(floor),high=Math.log10(Math.max(...positive,1)),span=Math.max(1,high-low),left=65,right=w-20,top=20,bottom=h-40;
 const px=i=>left+(right-left)*i/Math.max(1,history.length-1),py=v=>bottom-(Math.log10(Math.max(v,floor))-low)/span*(bottom-top);
 ctx.strokeStyle=theme().getPropertyValue('--line');ctx.fillStyle=theme().getPropertyValue('--muted');ctx.textAlign='right';
 for(let i=0;i<=4;i++){const y=top+(bottom-top)*i/4,value=10**(low+span*(1-i/4));ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText(value.toExponential(0),left-8,y+4);}
 ctx.textAlign='center';ctx.fillText('0',left,bottom+18);ctx.fillText(String(history.length-1),right,bottom+18);ctx.fillText('Solver iteration',w/2,h-6);
 ctx.strokeStyle='#ba8c30';ctx.setLineDash([5,4]);ctx.beginPath();ctx.moveTo(left,py(current.threshold));ctx.lineTo(right,py(current.threshold));ctx.stroke();ctx.setLineDash([]);
 ctx.strokeStyle=theme().getPropertyValue('--accent');ctx.lineWidth=2;ctx.beginPath();history.forEach((v,i)=>i?ctx.lineTo(px(i),py(v)):ctx.moveTo(px(i),py(v)));ctx.stroke();ctx.beginPath();ctx.arc(px(frame),py(history[frame]),5,0,Math.PI*2);ctx.fillStyle=ctx.strokeStyle;ctx.fill();ctx.lineWidth=1;
}
function updateFrame(value){if(!current)return;frame=Math.max(0,Math.min(current.iterates.length-1,value));$('frame').value=frame;$('step').textContent=`Iteration ${frame} / ${current.iterations}`;drawHeat();drawResidual();}
function setPlaying(value){playing=value;$('play').textContent=playing?'Pause':'Play';lastTime=0;}
function animate(time){if(playing&&current){if(!lastTime)lastTime=time;const duration=15000/Math.max(1,current.iterations);if(time-lastTime>=duration){const next=frame+Math.max(1,Math.floor((time-lastTime)/duration));lastTime=time;if(next>current.iterations){if($('loop').checked)updateFrame(0);else{updateFrame(current.iterations);setPlaying(false);}}else updateFrame(next);}}requestAnimationFrame(animate);}
requestAnimationFrame(animate);
$('play').onclick=()=>{if(!playing&&current&&frame===current.iterations)updateFrame(0);setPlaying(!playing);};
$('frame').oninput=()=>{setPlaying(false);updateFrame(Number($('frame').value));};
function config(){return {size:Number($('grid').value),contrast:Number($('contrast').value),limit:Number($('limit').value),jacobi:$('jacobi').checked};}
function fail(message){clearTimeout(watchdog);$('status').textContent=message;setPlaying(false);$('play').disabled=true;$('frame').disabled=true;worker?.terminate();worker=undefined;ready=false;if(workerURL){URL.revokeObjectURL(workerURL);workerURL=undefined;}}
function guard(){clearTimeout(watchdog);watchdog=setTimeout(()=>fail('The solver timed out. Change a setting to restart.'),30000);}
function dispatch(){if(!ready||!pending)return;guard();worker.postMessage(pending);pending=undefined;}
function start(){if(worker)return;if(!bundle.available){fail(bundle.reason);return;}workerURL=URL.createObjectURL(new Blob([bundle.worker_source],{type:'text/javascript'}));worker=new Worker(workerURL,{type:'module'});guard();worker.onerror=e=>fail('Solver failed: '+e.message);worker.onmessage=({data:m})=>{
 if(m.type==='ready'){clearTimeout(watchdog);ready=true;$('verification').textContent=`${m.verified} shared sparse fixtures passed in this browser. Live bundle SHA-256: ${bundle.sha256}`;dispatch();return;}
 if(m.type==='error'&&m.id===undefined){fail('Solver initialization failed: '+m.message);return;}
 if(m.id!==sequence)return;
 clearTimeout(watchdog);if(m.type==='error'){fail('Solver failed: '+m.message);return;}
 current=m.result;$('status').textContent=current.converged?'Converged: the true residual meets the requested tolerance.':`Stopped: ${current.reason.replaceAll('_',' ')}. The displayed estimate has not converged.`;
 $('stats').replaceChildren();for(const [value,label] of [[current.size**2,'unknowns'],[current.nnz,'nonzero coefficients'],[current.iterations,'accepted iterations'],[current.residuals.at(-1).toExponential(2),'final true residual']]){const e=add($('stats'),'div','');e.className='stat';add(e,'strong',value);add(e,'span',label);}
 $('frame').max=current.iterations;$('frame').disabled=false;$('play').disabled=current.iterations===0;setPlaying(false);updateFrame(current.iterations);
 };worker.postMessage({type:'init',bundle});}
function schedule(){sequence++;pending={type:'solve',id:sequence,config:config()};clearTimeout(timer);setPlaying(false);$('play').disabled=true;$('frame').disabled=true;$('status').textContent='Updating solver…';timer=setTimeout(()=>{start();dispatch();},300);}
for(const id of ['grid','contrast','limit','jacobi'])$(id).onchange=schedule;
function drawTimings(){const op=$('operation').value,rows=data.results.filter(r=>r.operation===op&&visible.has(r.implementation));$('timing-rows').replaceChildren();for(const r of [...rows].sort((a,b)=>a.implementation.localeCompare(b.implementation)||a.size-b.size)){const tr=add($('timing-rows'),'tr','');for(const v of [r.implementation,`${r.size} × ${r.size} / ${r.unknowns}`,r.status==='passed'?number(r.median_ns/1e6):r.status,r.status==='passed'?number(r.mad_ns/1e6):'—',r.cg_iterations??'—'])add(tr,'td',v);}
 const {ctx,w,h}=canvas('timings'),valid=rows.filter(r=>r.status==='passed'),maxX=Math.max(1,...valid.map(r=>r.unknowns)),maxY=Math.max(1e-9,...valid.map(r=>r.median_ns/1e6))*1.12,left=65,right=w-25,top=25,bottom=h-40;
 ctx.strokeStyle=theme().getPropertyValue('--line');ctx.fillStyle=theme().getPropertyValue('--muted');ctx.textAlign='right';for(let i=0;i<=4;i++){let y=bottom-(bottom-top)*i/4;ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText(number(maxY*i/4),left-8,y+4);}ctx.fillText('ms',left-8,14);ctx.textAlign='center';for(const n of [...new Set(valid.map(r=>r.unknowns))].sort((a,b)=>a-b))ctx.fillText(String(n),left+(right-left)*n/maxX,bottom+18);ctx.fillText('Unknowns',w/2,h-5);
 for(const [index,impl] of data.implementations.entries()){const points=valid.filter(r=>r.implementation===impl.id).sort((a,b)=>a.unknowns-b.unknowns);ctx.strokeStyle=ctx.fillStyle=colors[index%colors.length];ctx.lineWidth=2;ctx.beginPath();points.forEach((r,i)=>{const x=left+(right-left)*r.unknowns/maxX,y=bottom-(bottom-top)*r.median_ns/1e6/maxY;i?ctx.lineTo(x,y):ctx.moveTo(x,y);});ctx.stroke();for(const r of points){ctx.beginPath();ctx.arc(left+(right-left)*r.unknowns/maxX,bottom-(bottom-top)*r.median_ns/1e6/maxY,4,0,Math.PI*2);ctx.fill();}}
 if(!valid.length){ctx.fillStyle=theme().getPropertyValue('--muted');ctx.fillText('No included measurements',w/2,h/2);}
}
for(const [i,impl] of data.implementations.entries()){const label=add($('languages'),'label','');label.className='check';label.style.color=colors[i%colors.length];const input=document.createElement('input');input.type='checkbox';input.checked=true;input.onchange=()=>{input.checked?visible.add(impl.id):visible.delete(impl.id);drawTimings();};label.append(input,document.createTextNode(impl.name||impl.id));}
$('operation').onchange=drawTimings;
const storage=new Map(data.results.map(r=>[r.size,r]));for(const r of [...storage.values()].sort((a,b)=>a.size-b.size)){const tr=add($('storage'),'tr','');for(const v of [`${r.size} × ${r.size} / ${r.unknowns}`,number(r.nnz),number(r.logical_dense_bytes),number(r.logical_csr_bytes),number(r.logical_dense_bytes/r.logical_csr_bytes)+'×'])add(tr,'td',v);}
add($('method'),'p',`${data.created_at} · ${data.machine.os} ${data.machine.architecture} · ${data.machine.logical_cpus} logical CPUs`);
add($('method'),'p',`Git ${data.revision}${data.dirty?' (working tree contained changes)':''}`);add($('method'),'code',`Source SHA-256: ${data.source_sha256}`);
for(const [key,value] of Object.entries(data.methodology))add($('method'),'p',key.replaceAll('_',' ')+': '+value);
for(const impl of data.implementations){add($('toolchains'),'h3',`${impl.name||impl.id} · ${impl.status} · ${impl.checks?.filter(c=>c.passed).length||0} checks`);add($('toolchains'),'p',impl.toolchain||'');for(const command of impl.build_commands||[])add($('toolchains'),'code',command.join(' ')+'\n');}
function redraw(){drawHeat();drawResidual();drawTimings();}
let dark=matchMedia('(prefers-color-scheme: dark)').matches;try{dark=(localStorage.getItem('linear-a-theme')|| (dark?'dark':'light'))==='dark';}catch{}
function applyTheme(){document.documentElement.dataset.theme=dark?'dark':'light';$('theme').textContent=dark?'Light mode':'Dark mode';redraw();}
$('theme').onclick=()=>{dark=!dark;try{localStorage.setItem('linear-a-theme',dark?'dark':'light');}catch{}applyTheme();};
new ResizeObserver(redraw).observe(document.querySelector('main'));
window.addEventListener('pagehide',()=>{clearTimeout(timer);clearTimeout(watchdog);worker?.terminate();if(workerURL)URL.revokeObjectURL(workerURL);});
applyTheme();schedule();
