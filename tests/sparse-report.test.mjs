import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import {readFileSync} from 'node:fs';
const source=readFileSync(new URL('../benchmarks/sparse-report.mjs',import.meta.url),'utf8');
function setup(){
 const elements=new Map(),timers=new Map(),workers=[],animation=[],events={};let tid=0;
 const context=new Proxy({},{get:(o,k)=>o[k]??(()=>{})});
 const element=()=>({value:'',checked:false,disabled:false,textContent:'',children:[],style:{},append(...a){this.children.push(...a);},replaceChildren(){this.children=[];},getBoundingClientRect:()=>({width:500,height:340}),getContext:()=>context});
 const get=id=>{if(!elements.has(id))elements.set(id,element());return elements.get(id);};
 get('data').textContent=JSON.stringify({implementations:[],results:[],machine:{},methodology:{}});get('live').textContent=JSON.stringify({available:true,worker_source:'worker'});
 get('grid').value='16';get('contrast').value='2';get('limit').value='400';get('jacobi').checked=true;
 class Worker{constructor(){workers.push(this);this.sent=[];}postMessage(m){this.sent.push(m);}terminate(){this.terminated=true;}}
 const sandbox={document:{getElementById:get,createElement:element,createTextNode:t=>t,documentElement:{dataset:{}},querySelector:()=>element()},window:{addEventListener:(k,f)=>events[k]=f},Worker,Blob:class{},URL:{createObjectURL:()=> 'blob:test',revokeObjectURL(){}},devicePixelRatio:1,getComputedStyle:()=>({getPropertyValue:()=> '#000'}),ResizeObserver:class{observe(){}},matchMedia:()=>({matches:false}),localStorage:{getItem(){},setItem(){}},requestAnimationFrame:f=>animation.push(f),setTimeout:(f,delay)=>{timers.set(++tid,{f,delay});return tid;},clearTimeout:id=>timers.delete(id)};
 vm.runInNewContext(source,sandbox);
 const flush=delay=>{for(const [id,t] of [...timers])if(t.delay===delay){timers.delete(id);t.f();}};
 const ready=()=>{flush(300);workers.at(-1).onmessage({data:{type:'ready',verified:39}});};
 const result=id=>workers.at(-1).onmessage({data:{type:'result',id,result:{size:2,nnz:12,iterations:2,converged:true,reason:'converged',residuals:[10,1,.00001],threshold:.001,iterates:[[0,0,0,0],[1,2,3,4],[2,3,4,5]]}}});
 return {get,timers,workers,animation,events,flush,ready,result};
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
