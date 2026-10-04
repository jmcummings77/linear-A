import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import {readFileSync} from 'node:fs';
const source=readFileSync(new URL('../benchmarks/ic0-report.mjs',import.meta.url),'utf8');
function setup({width=900,height=320,recordedOverride}={}){
 const elements=new Map(),timers=new Map(),workers=[],frames=new Map(),events={},drawing=[];let id=0;
 const ctx=new Proxy({},{get:(o,k)=>o[k]??((...args)=>drawing.push([k,...args]))});
 const element=()=>({value:'',checked:false,disabled:false,textContent:'',children:[],width:700,height:700,style:{},handlers:{},getBoundingClientRect:()=>({width,height}),addEventListener(k,f){this.handlers[k]=f;},append(...v){this.children.push(...v);},replaceChildren(...v){this.children=v;},setAttribute(k,v){this[k]=v;},getContext:()=>ctx});
 const get=k=>{if(!elements.has(k))elements.set(k,element());return elements.get(k);};
 const recorded={implementations:[{id:'c',name:'C',status:'passed'}],created_at:'2026-10-03',revision:'abcdef',machine:{os:'test',architecture:'arm64'},methodology:{},results:[{implementation:'c',problem:'grid',size:12,operation:'ic0_total',ordering:'natural',status:'passed',median_ns:100,mad_ns:1},{implementation:'c',problem:'grid',size:12,operation:'cg',ordering:'natural',status:'failed',error:'zero pivot'}]};
 get('data').textContent=JSON.stringify(recordedOverride||recorded);get('live').textContent=JSON.stringify({available:true,capabilities:['ic0']});
 for(const [k,v] of Object.entries({language:'c',ordering:'natural','live-ordering':'natural','method-select':'ic0',problem:'grid',size:'12',metric:'total','live-size':'12','live-problem':'scrambled',contrast:'2',speed:'4',progress:'0'}))get(k).value=v;
 class Worker{constructor(){workers.push(this);this.sent=[];}postMessage(m){this.sent.push(m);}terminate(){this.terminated=true;}}
 const sandbox={document:{getElementById:get,createElement:element,createElementNS:element,body:{},documentElement:{dataset:{}}},window:{addEventListener:(k,f)=>events[k]=f},Worker,Blob:class{},URL:{createObjectURL:()=> 'blob:test',revokeObjectURL(){}},MutationObserver:class{observe(){}},devicePixelRatio:2,getComputedStyle:()=>({color:'#123',getPropertyValue:()=> '#ddd'}),performance:{now:()=>0},requestAnimationFrame:f=>{frames.set(++id,f);return id;},cancelAnimationFrame:i=>frames.delete(i),setTimeout:(f,delay)=>{timers.set(++id,{f,delay});return id;},clearTimeout:i=>timers.delete(i)};
 vm.runInNewContext(source,sandbox);
 const flush=delay=>{for(const [i,t] of [...timers])if(t.delay===delay){timers.delete(i);t.f();}};
 const reply=(id,failed=false,overrides={})=>{const r={ok:true,reason:'converged',iterations:2,setupMs:0,solveMs:1,nnz:3,frames:[[0,0,0,0],[1,1,1,1],[2,2,2,2]],residualFields:[[1,1,1,1],[.1,.1,.1,.1],[0,0,0,0]],residuals:[2,.2,0]};workers[0].onmessage({data:{type:'result',id,result:{n:4,size:2,ordering:'natural',orderingMs:0,runs:{cg:r,jacobi:r,ic0:failed?{ok:false,error:'nonpositive pivot',setupMs:0}:r,cholesky:r},...overrides}}});};
 return {get,workers,timers,frames,events,drawing,flush,reply};
}

test('IC0 controls debounce and ignore stale replies',()=>{
 const t=setup();assert.ok(t.get('saved-values').children.some(row=>row.children.some(cell=>String(cell.textContent).includes('zero pivot'))));
 t.workers[0].onmessage({data:{type:'ready'}});t.get('live-size').value='8';t.get('live-size').handlers.change();t.get('live-size').value='16';t.get('live-size').handlers.change();assert.equal([...t.timers.values()].filter(x=>x.delay===300).length,1);t.flush(300);assert.equal(t.workers[0].sent[1].config.size,16);
 t.reply(1);assert.equal(t.get('play').disabled,true);t.reply(3);assert.equal(t.get('play').disabled,false);assert.match(t.get('live-state').textContent,/Independent true-residual/);
});
test('IC0 animation loops, steps and clears timers on teardown',()=>{
 const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);
 t.get('next').handlers.click();assert.equal(t.get('progress').value,'1');assert.match(t.get('frame-code').textContent,/Accepted update 1/);
 t.get('loop').checked=true;t.get('play').handlers.click();let [id,f]=[...t.frames][0];t.frames.delete(id);f(11000);assert.equal(t.get('progress').value,'2');assert.equal(t.frames.size,1);
 t.get('progress').value='0';t.get('progress').handlers.input();assert.equal(t.frames.size,0);t.events.pagehide();assert.ok(t.workers[0].terminated);assert.equal(t.timers.size,0);
});
test('IC0 breakdown disables only unavailable playback; another method remains usable',()=>{
 const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1,true);assert.equal(t.get('play').disabled,true);assert.match(t.get('completed-status').textContent,/nonpositive pivot/);t.get('method-select').value='cg';t.get('method-select').handlers.change();assert.equal(t.get('play').disabled,false);
});
test('IC0 worker timeout and template IDs',()=>{const t=setup();t.flush(30000);assert.ok(t.workers[0].terminated);assert.match(t.get('live-state').textContent,/initialization timed out/);const html=readFileSync(new URL('../benchmarks/ic0-report.html',import.meta.url),'utf8'),ids=[...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]);assert.equal(new Set(ids).size,ids.length);});

test('changing the displayed method cannot enable stale playback while inputs are pending',()=>{const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);t.get('live-size').handlers.change();t.get('method-select').handlers.change();assert.equal(t.get('play').disabled,true);});

test('saved chart uses display-width coordinates and readable labels on mobile',()=>{
 const t=setup({width:300}),chart=t.get('saved-chart');
 const [,,width,height]=chart.viewBox.split(' ').map(Number);
 assert.equal(width,300);assert.equal(chart.style.height,height+'px');
 const texts=chart.children.filter(node=>node['font-size']!==undefined);
 assert.ok(texts.length>=8);assert.ok(texts.every(node=>node['font-size']===13));
 const label=texts.find(node=>node.textContent==='IC(0)-CG total'),value=texts.find(node=>node.textContent==='100 ns');
 assert.equal(label.y,value.y);assert.equal(value['text-anchor'],'end');assert.equal(value.x,300);
 const bars=chart.children.filter(node=>node.rx===3);assert.equal(Math.max(...bars.map(node=>node.width)),300);
 assert.ok(t.get('saved-values').children.some(row=>row.children.some(cell=>String(cell.textContent).includes('zero pivot'))));
});

function capturedRun(residuals,overrides={}){
 return {ok:true,reason:'converged',iterations:residuals.length-1,setupMs:0,solveMs:1,nnz:3,
  frames:residuals.map((_,index)=>Array(4).fill(index)),residualFields:residuals.map(value=>Array(4).fill(value/2)),residuals,...overrides};
}

test('convergence fits growth and very small residuals, and labels exact zeros separately',()=>{
 const t=setup({width:300});t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.drawing.length=0;
 t.reply(1,false,{runs:{cg:capturedRun([2,20,1e-16,0]),jacobi:capturedRun([2,2e-3,0]),ic0:capturedRun([2,1e-8,0]),cholesky:capturedRun([2,0],{iterations:0})}});
 const labels=t.drawing.filter(([operation])=>operation==='fillText').map(([,text])=>text);
 assert.ok(labels.includes('1e2'));assert.ok(labels.includes('1e-16'));assert.ok(labels.includes('zero'));
 assert.equal(t.get('convergence').width,600);assert.equal(t.get('convergence').style.height,'300px');
 assert.ok(t.drawing.some(([operation,...dash])=>operation==='setLineDash'&&String(dash[0])==='9,5'));
 assert.match(t.get('convergence-note').textContent,/2\.000e-9/);assert.match(t.get('convergence-note').textContent,/separate zero row/);
 for(const [operation,...values] of t.drawing)if(['moveTo','lineTo','arc','rect'].includes(operation))assert.ok(values.every(Number.isFinite),`${operation} must have finite coordinates`);
 // The CG residual grows from 2 to 20; its second circular marker must move up.
 const circles=t.drawing.filter(([operation])=>operation==='arc');assert.ok(circles[1][2]<circles[0][2]);
});

test('current-frame residual is separate from completed outcome and field scales stay common',()=>{
 const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);
 t.reply(1,false,{runs:{cg:capturedRun([2,20,0],{frames:[[0,0,0,0],[50,0,0,0],[2,2,2,2]]}),jacobi:capturedRun([2,0]),ic0:capturedRun([2,.2,0]),cholesky:capturedRun([2,0],{iterations:0})}});
 assert.equal(t.get('current-residual').textContent,'2.0000e+0');assert.equal(t.get('completed-status').textContent,'converged after 2 updates');
 assert.match(t.get('frame-note').textContent,/Accepted update 0/);assert.doesNotMatch(t.get('frame-note').textContent,/converged/);
 assert.equal(t.get('solution-max').textContent,'+50.0');assert.equal(t.get('residual-max').textContent,'+10.0');
 t.get('next').handlers.click();assert.equal(t.get('current-residual').textContent,'2.0000e-1');assert.equal(t.get('current-threshold').textContent,'2.0000e-9');
 t.get('method-select').value='cg';t.get('method-select').handlers.change();
 assert.equal(t.get('solution-max').textContent,'+50.0');assert.equal(t.get('residual-max').textContent,'+10.0');
 const html=readFileSync(new URL('../benchmarks/ic0-report.html',import.meta.url),'utf8');
 assert.match(html,/<p id="frame-note" aria-live="off">/);assert.match(html,/<dl class="ic0-current" aria-live="off">/);
});
