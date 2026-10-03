import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import {readFileSync} from 'node:fs';
const source=readFileSync(new URL('../benchmarks/cholesky-report.mjs',import.meta.url),'utf8');
function setup(){
 const elements=new Map(),timers=new Map(),workers=[],frames=new Map(),events={};let id=0;
 const ctx=new Proxy({},{get:(o,k)=>o[k]??(()=>{})});
 const element=()=>({value:'',checked:false,disabled:false,textContent:'',children:[],width:700,height:700,handlers:{},addEventListener(k,f){this.handlers[k]=f;},append(...v){this.children.push(...v);},replaceChildren(...v){this.children=v;},setAttribute(k,v){this[k]=v;},getContext:()=>ctx});
 const get=k=>{if(!elements.has(k))elements.set(k,element());return elements.get(k);};
 const recorded={implementations:[{id:'c',name:'C',status:'passed'}],created_at:'2026-10-03',revision:'abcdef',machine:{os:'test',architecture:'arm64'},methodology:{},results:[{implementation:'c',problem:'grid',size:12,operation:'total',ordering:'natural',status:'passed',median_ns:100,mad_ns:1},{implementation:'c',problem:'grid',size:12,operation:'total',ordering:'rcm',status:'failed',error:'zero pivot'}]};
 get('data').textContent=JSON.stringify(recorded);get('live').textContent=JSON.stringify({available:true,capabilities:['cholesky']});
 for(const [k,v] of Object.entries({language:'c',problem:'grid',size:'12',metric:'total','live-size':'12','live-problem':'scrambled',contrast:'2',speed:'4',progress:'0'}))get(k).value=v;
 class Worker{constructor(){workers.push(this);this.sent=[];}postMessage(m){this.sent.push(m);}terminate(){this.terminated=true;}}
 const sandbox={document:{getElementById:get,createElement:element,createElementNS:element,body:{},documentElement:{}},window:{addEventListener:(k,f)=>events[k]=f},Worker,Blob:class{},URL:{createObjectURL:()=> 'blob:test',revokeObjectURL(){}},MutationObserver:class{observe(){}},getComputedStyle:()=>({color:'#123'}),performance:{now:()=>0},requestAnimationFrame:f=>{frames.set(++id,f);return id;},cancelAnimationFrame:i=>frames.delete(i),setTimeout:(f,delay)=>{timers.set(++id,{f,delay});return id;},clearTimeout:i=>timers.delete(i)};
 vm.runInNewContext(source,sandbox);
 const flush=delay=>{for(const [i,t] of [...timers])if(t.delay===delay){timers.delete(i);t.f();}};
 const reply=(id)=>{const solver={ok:true,points:[[0,0,-1,2],[1,0,0,1],[1,1,-1,2]],nnz:3,fillCount:1,residual:1e-14,reconstructionError:1e-14,symbolicMilliseconds:0,factorMilliseconds:1,solveMilliseconds:1,totalMilliseconds:2};workers[0].onmessage({data:{type:'result',id,result:{n:2,points:[[0,0,1],[1,1,2]],p:[1,0],inverse:[1,0],width:0,newWidth:0,natural:solver,rcm:solver}}});};
 return {get,workers,timers,frames,events,flush,reply};
}
test('Cholesky controls debounce, ignore stale replies and retain recorded failures',()=>{
 const t=setup();assert.ok(t.get('saved-values').children.some(row=>row.children.some(cell=>String(cell.textContent).includes('zero pivot'))));
 t.workers[0].onmessage({data:{type:'ready'}});t.get('live-size').value='8';t.get('live-size').handlers.change();t.get('live-size').value='16';t.get('live-size').handlers.change();
 assert.equal([...t.timers.values()].filter(x=>x.delay===300).length,1);t.flush(300);assert.equal(t.workers[0].sent[1].config.size,16);
 t.reply(1);assert.equal(t.get('play').disabled,true);t.reply(3);assert.equal(t.get('play').disabled,false);assert.match(t.get('live-state').textContent,/Both factorizations finished/);
});
test('animation loops, scrubbing stops playback, teardown clears timers and workers',()=>{
 const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);
 t.get('loop').checked=true;t.get('play').handlers.click();let [id,f]=[...t.frames][0];t.frames.delete(id);f(11000);assert.equal(t.get('progress').value,'1000');assert.equal(t.frames.size,1);
 t.get('progress').value='300';t.get('progress').handlers.input();assert.equal(t.frames.size,0);assert.equal(t.get('play').textContent,'Animate elimination');
 t.events.pagehide();assert.ok(t.workers[0].terminated);assert.equal(t.timers.size,0);
});
test('worker timeout is visible and terminates work',()=>{const t=setup();t.flush(30000);assert.ok(t.workers[0].terminated);assert.match(t.get('live-state').textContent,/initialization timed out/);});
test('HTML IDs are unique, including executable data and live section',()=>{
 const html=readFileSync(new URL('../benchmarks/cholesky-report.html',import.meta.url),'utf8'),ids=[...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]);assert.equal(new Set(ids).size,ids.length);
});

test('pivot buttons step the actual elimination history',()=>{const t=setup();t.workers[0].onmessage({data:{type:'ready'}});t.flush(300);t.reply(1);t.get('next').handlers.click();assert.equal(t.get('progress').value,'500');assert.match(t.get('frame-code').textContent,/pivot k = 0/);t.get('previous').handlers.click();assert.equal(t.get('progress').value,'0');});
