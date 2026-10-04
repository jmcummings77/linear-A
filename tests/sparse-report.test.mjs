import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import {readFileSync} from 'node:fs';
const source=readFileSync(new URL('../benchmarks/sparse-report.mjs',import.meta.url),'utf8');
function setup(recorded={implementations:[],results:[],machine:{},methodology:{}},live={available:true,worker_source:'worker',capabilities:['gmres','ilu0','rhs_cache']},dimensions={width:500,height:340}){
 const elements=new Map(),timers=new Map(),workers=[],animation=[],events={},drawing=[];let tid=0;
 const context=new Proxy({},{get:(o,k)=>o[k]??((...args)=>drawing.push([k,...args]))});
 const element=()=>({value:'',checked:false,disabled:false,textContent:'',children:[],style:{},setAttribute(name,value){this[name]=value;},removeAttribute(name){delete this[name];},append(...a){this.children.push(...a);},replaceChildren(){this.children=[];},getBoundingClientRect:()=>dimensions,getContext:()=>context});
 const get=id=>{if(!elements.has(id))elements.set(id,element());return elements.get(id);};
 get('data').textContent=JSON.stringify(recorded);get('live').textContent=JSON.stringify(live);
 get('operation').value='csr_spmv';get('timing-scale').value='linear';
 get('solver').value='gmres';get('restart').value='20';get('diffusivity').value='.2';get('speed').value='4';get('angle').value='30';
 get('grid').value='16';get('contrast').value='2';get('limit').value='400';get('preconditioner').value='jacobi';get('hot').value='100';get('rhs-count').value='8';
 class Worker{constructor(){workers.push(this);this.sent=[];}postMessage(m){this.sent.push(m);}terminate(){this.terminated=true;}}
 const sandbox={document:{getElementById:get,createElement:element,createTextNode:t=>t,documentElement:{dataset:{}},querySelector:()=>element()},window:{addEventListener:(k,f)=>events[k]=f},Worker,Blob:class{},URL:{createObjectURL:()=> 'blob:test',revokeObjectURL(){}},devicePixelRatio:1,getComputedStyle:()=>({getPropertyValue:()=> '#000'}),ResizeObserver:class{observe(){}},matchMedia:()=>({matches:false}),localStorage:{getItem(){},setItem(){}},requestAnimationFrame:f=>animation.push(f),setTimeout:(f,delay)=>{timers.set(++tid,{f,delay});return tid;},clearTimeout:id=>timers.delete(id)};
 vm.runInNewContext(source,sandbox);
 const flush=delay=>{for(const [id,t] of [...timers])if(t.delay===delay){timers.delete(id);t.f();}};
 const ready=()=>{flush(300);workers.at(-1).onmessage({data:{type:'ready',verified:39}});};
 const result=(id,overrides={})=>workers.at(-1).onmessage({data:{type:'result',id,result:{size:2,nnz:12,iterations:2,converged:true,reason:'converged',residuals:[10,1,.00001],threshold:.001,iterates:[[0,0,0,0],[1,2,3,4],[2,3,4,5]],...overrides}}});
 return {get,timers,workers,animation,events,drawing,flush,ready,result};
}
test('debounce coalesces rapid changes and ignores stale worker results',()=>{
 const t=setup();t.get('grid').value='8';t.get('grid').onchange();t.get('grid').value='24';t.get('grid').onchange();
 assert.equal([...t.timers.values()].filter(x=>x.delay===300).length,1);t.ready();assert.equal(t.workers[0].sent[1].config.size,24);
 t.result(1);assert.equal(t.get('play').disabled,true);t.result(3);assert.equal(t.get('play').disabled,false);assert.match(t.get('status').textContent,/Converged/);
});
test('playback advances captured iterations, loops and pauses for manual stepping',()=>{
 const t=setup();t.ready();t.result(1);t.get('play').onclick();assert.equal(t.get('play').textContent,'Pause');assert.equal(t.get('frame').value,0);
 t.animation.shift()(100);t.animation.shift()(8000);assert.equal(t.get('frame').value,1);
 t.get('loop').checked=true;t.animation.shift()(24000);assert.equal(t.get('frame').value,0);
 t.get('frame').value='2';t.get('frame').oninput();assert.equal(t.get('play').textContent,'Play');assert.equal(t.get('frame').value,2);
});
test('initialization errors and watchdogs terminate workers and allow restart',()=>{
 const t=setup();t.flush(300);t.workers[0].onmessage({data:{type:'error',message:'bad bundle'}});assert.match(t.get('status').textContent,/initialization failed: bad bundle/);assert.ok(t.workers[0].terminated);
 t.get('grid').onchange();t.flush(300);assert.equal(t.workers.length,2);t.flush(30000);assert.ok(t.workers[1].terminated);assert.match(t.get('status').textContent,/timed out/);
 t.get('grid').onchange();t.flush(300);t.events.pagehide();assert.ok(t.workers[2].terminated);assert.equal(t.timers.size,0);
});

const recorded={
 implementations:[{id:'csharp',name:'C#',status:'passed'},{id:'python',name:'Python',status:'passed'}],
 machine:{os:'macOS',architecture:'arm64'},methodology:{timing:'Output allocation included.'},created_at:'2026-08-05T12:00:00Z',revision:'1234567890abcdef',
 results:[
  {implementation:'csharp',operation:'csr_spmv',size:8,unknowns:64,nnz:288,logical_dense_bytes:32768,logical_csr_bytes:5128,status:'passed',median_ns:2000000,mad_ns:1000},
  {implementation:'csharp',operation:'csr_spmv',size:16,unknowns:256,nnz:1216,logical_dense_bytes:524288,logical_csr_bytes:21512,status:'passed',median_ns:4000000,mad_ns:8000},
  {implementation:'python',operation:'csr_spmv',size:8,unknowns:64,nnz:288,logical_dense_bytes:32768,logical_csr_bytes:5128,status:'passed',median_ns:80000000,mad_ns:2000},
  {implementation:'python',operation:'csr_spmv',size:16,unknowns:256,nnz:1216,logical_dense_bytes:524288,logical_csr_bytes:21512,status:'failed'}
 ]
};
test('saved timing controls preserve values and expose failed rows while changing chart scale',()=>{
 const original=JSON.stringify(recorded),t=setup(recorded);
 const rows=()=>t.get('timing-rows').children.map(row=>row.children.map(cell=>cell.textContent));
 assert.deepEqual(rows()[0],['C#','8 × 8 / 64','2 ms','1 µs','—','—','—']);
 assert.equal(rows()[3][2],'failed');
 assert.ok(t.drawing.some(([op,text])=>op==='fillText'&&text==='80'));
 t.drawing.length=0;t.get('timing-scale').value='log';t.get('timing-scale').onchange();
 assert.match(t.get('timing-note').textContent,/Logarithmic/);
 assert.deepEqual(rows()[0],['C#','8 × 8 / 64','2 ms','1 µs','—','—','—']);
 assert.ok(t.drawing.some(([op,text])=>op==='fillText'&&text==='100'));
 assert.equal(JSON.stringify(recorded),original);
});
test('visible legend and language highlight agree when a hidden series is selected',()=>{
 const t=setup(recorded),labels=t.get('languages').children;
 assert.equal(labels.length,2);assert.equal(labels[0].children[2],'C#');
 const pythonInput=labels[1].children[0];pythonInput.checked=false;pythonInput.onchange();
 assert.equal(t.get('timing-rows').children.length,2);assert.match(t.get('timing-note').textContent,/1\/2 languages/);
 t.get('highlight').value='python';t.get('highlight').onchange();
 assert.equal(pythonInput.checked,true);assert.equal(t.get('timing-rows').children.length,4);
 assert.equal(t.get('timing-rows').children[2].className,'series-selected');
 assert.match(t.get('timing-note').textContent,/Highlighting Python/);
 pythonInput.checked=false;pythonInput.onchange();assert.equal(t.get('highlight').value,'');
});
test('logarithmic timing view keeps zero measurements in its table and reports omitted points',()=>{
 const data=structuredClone(recorded);data.results[0].median_ns=0;
 const t=setup(data);t.get('timing-scale').value='log';t.get('timing-scale').onchange();
 assert.equal(t.get('timing-rows').children[0].children[2].textContent,'0 ns');
 assert.match(t.get('timing-note').textContent,/1 zero timing/);
 for(const [operation,...args] of t.drawing)if(['moveTo','lineTo','arc','rect'].includes(operation))assert.ok(args.every(Number.isFinite));
 t.get('languages').children.forEach(label=>{const input=label.children[0];input.checked=false;input.onchange();});
 assert.ok(t.drawing.some(([op,text])=>op==='fillText'&&text==='No included measurements'));
 assert.equal(typeof t.events['matrix-theme-change'],'function');
});

test('CG disables advection and GMRES changes debounce together',()=>{
 const t=setup();t.get('solver').value='cg';t.get('solver').onchange();
 assert.equal(t.get('speed').value,'0');assert.equal(t.get('restart').disabled,true);assert.equal(t.get('angle').disabled,true);
 t.get('solver').value='gmres';t.get('solver').onchange();t.get('restart').value='5';t.get('restart').onchange();
 t.get('angle').value='90';t.get('angle').oninput();t.ready();
 const config=t.workers[0].sent.at(-1).config;assert.equal(config.solver,'gmres');assert.equal(config.restart,5);assert.equal(config.angle,90);
 assert.equal(t.get('restart').disabled,false);
});
test('GMRES timing rows include accepted work and logical workspace',()=>{
 const data=structuredClone(recorded);data.results=data.results.slice(0,1);Object.assign(data.results[0],{solver_iterations:12,logical_workspace_bytes:4096});
 const t=setup(data),cells=t.get('timing-rows').children[0].children;
 assert.equal(cells[4].textContent,12);assert.equal(cells[5].textContent,'4 KiB');
});

test('ILU selection and RHS changes debounce; CG removes incompatible ILU',()=>{
 const t=setup();t.get('preconditioner').value='ilu0';t.get('preconditioner').onchange();t.get('hot').value='50';t.get('hot').oninput();t.ready();
 let config=t.workers[0].sent.at(-1).config;assert.equal(config.ilu,true);assert.equal(config.jacobi,false);assert.equal(config.hot,50);
 t.get('solver').value='cg';t.get('solver').onchange();assert.equal(t.get('preconditioner').value,'jacobi');assert.equal(t.get('ilu-option').disabled,true);
});
test('reuse cost model includes setup once and exposes the crossover estimate',()=>{
 const data=structuredClone(recorded);data.suite='ilu-reuse-v1';data.implementations=data.implementations.slice(0,1);
 data.results=Object.entries({ilu_setup:100000,gmres_ilu_reused:10000,gmres_ilu_total:110000,gmres_none:100000,gmres_jacobi:30000}).map(([operation,median_ns])=>({...recorded.results[0],operation,median_ns}));
 const t=setup(data);assert.equal(t.get('reuse-section').hidden,false);assert.match(t.get('reuse-note').textContent,/after 6 RHS/);
 t.get('rhs-count').value='16';t.get('rhs-count').oninput();assert.match(t.get('reuse-note').textContent,/setup \+ 16/);
 assert.ok(t.drawing.some(([op,text])=>op==='fillText'&&text==='260 µs'));
});

test('preserved legacy live bundles do not advertise unsupported solver controls',()=>{
 const t=setup(undefined,{available:true,worker_source:'old worker',fixtures:[{op:'cg'}]});
 assert.equal(t.get('solver').value,'cg');assert.equal(t.get('solver').disabled,true);
 assert.equal(t.get('ilu-option').disabled,true);assert.equal(t.get('hot').disabled,true);assert.equal(t.get('capability-note').hidden,false);
 t.ready();assert.equal(t.workers[0].sent.at(-1).config.speed,0);
});


test('guided presets choose supported controls and preserve legacy capability restrictions',()=>{
 const t=setup();t.get('live-preset').value='ilu';t.get('live-preset').onchange();t.ready();
 let config=t.workers[0].sent.at(-1).config;assert.equal(config.solver,'gmres');assert.equal(config.ilu,true);assert.equal(config.speed,4);
 t.get('live-preset').value='diffusion';t.get('live-preset').onchange();t.flush(300);
 config=t.workers[0].sent.at(-1).config;assert.equal(config.solver,'cg');assert.equal(config.speed,0);assert.equal(config.ilu,false);assert.equal(config.jacobi,true);
 assert.equal(t.get('restart').disabled,true);
 t.get('live-preset').value='transport';t.get('live-preset').onchange();t.flush(300);
 config=t.workers[0].sent.at(-1).config;assert.equal(config.solver,'gmres');assert.equal(config.speed,4);assert.equal(t.get('restart').disabled,false);
 t.get('grid').value='24';t.get('grid').onchange();assert.equal(t.get('live-preset').value,'custom');
 const legacy=setup(undefined,{available:true,worker_source:'old worker',fixtures:[{op:'cg'}]});
 assert.equal(legacy.get('preset-transport').disabled,true);assert.equal(legacy.get('preset-ilu').disabled,true);
 legacy.get('live-preset').value='transport';legacy.get('live-preset').onchange();legacy.ready();
 assert.equal(legacy.workers[0].sent.at(-1).config.solver,'cg');assert.equal(legacy.get('hot').disabled,true);
 const gmres=setup(undefined,{available:true,worker_source:'gmres worker',capabilities:['gmres']});
 assert.equal(gmres.get('preset-transport').disabled,false);assert.equal(gmres.get('preset-ilu').disabled,true);
 gmres.get('live-preset').value='ilu';gmres.get('live-preset').onchange();gmres.ready();assert.equal(gmres.workers[0].sent.at(-1).config.ilu,false);
});

test('scrubbing, playback and table selection expose the same current residual and restart position',()=>{
 const t=setup();t.ready();t.result(1,{solver:'gmres',estimatedResiduals:[10,.8,.000009],restarts:[1]});
 assert.equal(t.get('frame-residual').textContent,'1.000e-5');assert.equal(t.get('frame-estimate').textContent,'9.000e-6');assert.equal(t.get('frame-tolerance').textContent,'1.000e-3');
 assert.equal(t.get('frame-restart').textContent,'1 step after restart 1');
 const rows=t.get('iteration-rows').children;assert.equal(rows.length,3);assert.equal(rows[2]['aria-current'],'step');
 t.get('frame').value='1';t.get('frame').oninput();assert.equal(t.get('frame-residual').textContent,'1.000e+0');assert.equal(t.get('frame-estimate').textContent,'8.000e-1');assert.equal(t.get('frame-restart').textContent,'Restart at iteration 1');assert.equal(rows[2]['aria-current'],undefined);assert.equal(rows[1]['aria-current'],'step');
 rows[0].children[0].children[0].onclick();assert.equal(t.get('frame').value,0);assert.equal(t.get('play').textContent,'Play');assert.equal(t.get('frame-restart').textContent,'Initial guess');
 t.get('play').onclick();t.animation.shift()(100);t.animation.shift()(8000);assert.equal(t.get('frame').value,1);assert.equal(t.get('frame-residual').textContent,rows[1].children[1].textContent);
 const cg=setup();cg.ready();cg.result(1);assert.equal(cg.get('frame-estimate').textContent,'Not applicable · CG');assert.equal(cg.get('frame-restart').textContent,'Not used by CG');
});

test('reuse chart keeps the full plotting width on a narrow screen',()=>{
 const data=structuredClone(recorded);data.suite='ilu-reuse-v1';data.implementations=data.implementations.slice(0,1);
 data.results=Object.entries({ilu_setup:100000,gmres_ilu_reused:10000,gmres_ilu_total:110000,gmres_none:100000,gmres_jacobi:30000}).map(([operation,median_ns])=>({...recorded.results[0],operation,median_ns}));
 const t=setup(data,undefined,{width:280,height:300});
 assert.equal(t.get('reuse-nav').hidden,false);
 const bars=t.drawing.filter(([op])=>op==='fillRect');assert.equal(bars.length,4);
 assert.ok(bars.every(([,x,y,width])=>x===16&&y>=34&&width>0&&width<=248));
 assert.equal(Math.max(...bars.map(([,x,y,width])=>width)),248);
 const totals=t.get('reuse-values').children.map(row=>row.textContent);assert.ok(totals.includes('ILU reused: 180 µs'));
});
