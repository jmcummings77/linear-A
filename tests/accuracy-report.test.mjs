import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';
import {test} from 'node:test';
const source=readFileSync(new URL('../benchmarks/accuracy-report.mjs',import.meta.url),'utf8');
function harness(available=true){
 class Element {
  constructor(id='',tag='DIV'){this.id=id;this.tagName=tag;this.value='';this.textContent='';this.children=[];this.events={};this.clientWidth=640;this.attrs={};}
  get value(){return this._value;}
  set value(value){this._value=String(value);}
  addEventListener(type,fn){(this.events[type]||=[]).push(fn);}
  dispatch(type){for(const fn of this.events[type]||[])fn();}
  append(child){this.children.push(child);}
  setAttribute(key,value){this.attrs[key]=value;}
  getContext(){return new Proxy({}, {get:()=>()=>{},set:()=>true});}
 }
 const nodes=new Map(),get=id=>{if(!nodes.has(id))nodes.set(id,new Element(id));return nodes.get(id);};
 for(const [id,tag,value] of [['preset','SELECT','stable'],['algorithm','SELECT','lu'],['target','SELECT','5'],['delta','INPUT','0.0000001']])Object.assign(get('accuracy-'+id),{tagName:tag,value});
 get('accuracy-panel').querySelectorAll=()=>['preset','algorithm','target','delta'].map(id=>get('accuracy-'+id));
 get('live-data').textContent=JSON.stringify({available,accuracy_worker_source:'export {};'});
 const timers=new Map(),workers=[],revoked=[];let clock=0,next=0;
 class Worker {constructor(){this.sent=[];workers.push(this);}postMessage(message){this.sent.push(message);}terminate(){this.terminated=true;}}
 const sandbox={document:{getElementById:get,createElement:tag=>new Element('',tag.toUpperCase()),documentElement:{}},Worker,WebAssembly:{},Blob:class{},URL:{createObjectURL:()=> 'blob:test',revokeObjectURL:value=>revoked.push(value)},MutationObserver:class{observe(){}},devicePixelRatio:1,getComputedStyle:()=>({getPropertyValue:()=> '#888'}),addEventListener(){},setTimeout(fn,ms){const id=++next;timers.set(id,{fn,due:clock+ms});return id;},clearTimeout:id=>timers.delete(id)};
 vm.runInNewContext(source,sandbox);
 function advance(ms){clock+=ms;for(const [id,timer] of [...timers])if(timer.due<=clock){timers.delete(id);timer.fn();}}
 return {get,workers,advance,revoked};
}
const response=id=>({data:{id,type:'accuracy',checks:{passed:166,total:166},result:{a:[4,1,1,3],b:[6,7],perturbed_a:[4,1,1,3],perturbed_b:[6,7],base:{x:[1,2],reciprocal_condition:.44,residual_infinity:0,backward_error:0},perturbed:{x:[1,2],reciprocal_condition:.44,residual_infinity:0,backward_error:0},actual_delta:0,relative_solution_change:0}}});
test('rapid edits debounce, stale answers are ignored, and worker is reused',()=>{
 const h=harness();h.advance(299);assert.equal(h.workers.length,0);h.advance(1);
 const worker=h.workers[0],first=worker.sent[0];
 worker.onmessage(response(first.id));assert.match(h.get('accuracy-solution').textContent,/1.00000/);
 const input=h.get('accuracy-a').children[0];input.value='5';input.dispatch('input');h.advance(200);input.value='6';input.dispatch('input');
 worker.onmessage(response(first.id));assert.equal(h.get('accuracy-solution').textContent,'—');
 h.advance(299);assert.equal(worker.sent.length,1);h.advance(1);
 assert.equal(h.workers.length,1);assert.equal(worker.sent.length,2);assert.equal(worker.sent[1].config.a[0],6);
 worker.onmessage(response(worker.sent[1].id));assert.match(h.get('accuracy-status').textContent,/166\/166/);
});
test('invalid edits clear old output and missing WASM never starts a worker',()=>{
 const h=harness();h.advance(300);const input=h.get('accuracy-b').children[0];input.value='';input.dispatch('input');h.advance(300);
 assert.equal(h.workers[0].sent.length,1);assert.match(h.get('accuracy-status').textContent,/finite matrix entries/);
 const unavailable=harness(false);unavailable.advance(300);assert.equal(unavailable.workers.length,0);assert.match(unavailable.get('accuracy-status').textContent,/Regenerate/);
});
test('timeout terminates worker, clears results, and editing can retry',()=>{
 const h=harness();h.advance(300);h.advance(20000);assert.equal(h.workers[0].terminated,true);assert.match(h.get('accuracy-status').textContent,/timed out/);
 h.get('accuracy-delta').dispatch('input');h.advance(300);assert.equal(h.workers.length,2);
});
