import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import vm from 'node:vm';

const source=readFileSync(new URL('../../benchmarks/live-report.mjs',import.meta.url),'utf8');
const template=readFileSync(new URL('../../benchmarks/report.html',import.meta.url),'utf8');
const controlIds=['live-check','live-run','live-operation','live-size','live-samples'];

class Element {
  constructor(tag){this.tagName=tag.toUpperCase();this.children=[];this.parentElement=null;this.attributes={};this._value='';this._text='';this.selectedIndex=-1;this.disabled=false;this.hidden=false;this.open=false;}
  set value(value){if(this.tagName==='SELECT')this.selectedIndex=this.options.findIndex(option=>option.value===String(value));else this._value=String(value);}
  get value(){return this.tagName==='SELECT'?(this.options[this.selectedIndex]?.value??''):this._value;}
  get options(){return this.children.filter(child=>child.tagName==='OPTION');}
  set textContent(value){this.children=[];this._text=String(value);}
  get textContent(){return this._text+this.children.map(child=>child.textContent).join('');}
  append(...children){for(const child of children){child.parentElement=this;this.children.push(child);if(this.tagName==='SELECT'&&this.selectedIndex<0)this.selectedIndex=0;}}
  replaceChildren(...children){this.children=[];this._text='';this.selectedIndex=-1;this.append(...children);}
  setAttribute(name,value){this.attributes[name]=String(value);}
  click(){if(!this.disabled)this.onclick?.();}
  remove(){if(this.parentElement)this.parentElement.children=this.parentElement.children.filter(child=>child!==this);this.parentElement=null;}
}

function page(){
  const nodes=new Map(),workers=[],created=[],revoked=[],timers=new Map();let sequence=0;
  const attributes=raw=>Object.fromEntries([...raw.matchAll(/([\w-]+)(?:="([^"]*)")?/g)].map(match=>[match[1],match[2]??'']));
  for(const match of template.matchAll(/<([a-z][\w-]*)\b([^>]*\bid="([^"]+)"[^>]*)>/g)){
    const [,tag,raw,id]=match,attrs=attributes(raw),node=new Element(tag);nodes.set(id,node);
    node.disabled='disabled'in attrs;node.hidden='hidden'in attrs;node.attributes=attrs;
    if(tag==='select'){
      const inner=template.slice(match.index+match[0].length).split('</select>')[0];
      for(const optionMatch of inner.matchAll(/<option\b([^>]*)>([^<]*)<\/option>/g)){
        const option=new Element('option'),optionAttrs=attributes(optionMatch[1]);option.value=optionAttrs.value??optionMatch[2];option.textContent=optionMatch[2];node.append(option);if('selected'in optionAttrs)node.selectedIndex=node.options.length-1;
      }
    }
  }
  const get=id=>{assert.ok(nodes.has(id),`unknown DOM id ${id}`);return nodes.get(id);};
  const fixtures=[{name:'rotation3d axis-angle'},{name:'rectangular multiplication'}];
  get('live-data').textContent=JSON.stringify({available:true,sha256:'ab'.repeat(32),byte_length:80,fixtures,worker_source:'export {};'});
  class Worker {
    constructor(url,options){this.url=url;this.options=options;this.messages=[];this.terminations=0;workers.push(this);}
    postMessage(message){this.messages.push(structuredClone(message));}
    terminate(){this.terminations++;}
    deliver(message){this.onmessage({data:structuredClone(message)});}
    fail(message){this.onerror({message,preventDefault(){}});}
  }
  const body=new Element('body');
  vm.runInNewContext(source,{document:{getElementById:get,createElement:tag=>new Element(tag),body},navigator:{userAgent:'Mozilla/5.0 Firefox/145.0'},Worker,WebAssembly:{},Blob,URL:{createObjectURL(blob){const url=`blob:live-${++sequence}`;created.push({url,blob});return url;},revokeObjectURL:url=>revoked.push(url)},setTimeout(callback,delay){const id=++sequence;timers.set(id,{callback,delay});return id;},clearTimeout:id=>timers.delete(id)});
  return {get,workers,created,revoked,timers,
    choose(id,value){assert.equal(get(id).disabled,false,`${id} must be usable`);get(id).value=value;assert.equal(get(id).value,String(value),`${value} must be a real select option`);get(id).onchange?.();},
    click(id){get(id).click();},
    export(){assert.equal(get('live-export').disabled,false);get('live-export').click();return JSON.parse(get('live-export-json').textContent);},
    snapshot(){return JSON.stringify({status:get('live-status').textContent,statusClass:get('live-status').className,summary:get('live-check-summary').textContent,checks:get('live-check-rows').textContent,results:get('live-result-rows').textContent,disabled:controlIds.map(id=>get(id).disabled),busy:get('live-panel').attributes['aria-busy'],stop:get('live-stop').disabled,export:get('live-export').disabled});},
  };
}

function checks(){return {passed:2,total:2,checks:[{name:'rotation3d axis-angle',passed:true},{name:'rectangular multiplication',passed:true}]};}
function benchmark(samples,checksumVerified=true){
  return {operation:'rotation3d',size:3,seed:17,iterations:7,samples:Array.from({length:samples},(_,index)=>({elapsed_ms:.14+index*.007,ns_per_op:20000+index*1000,checksum:11.75})),median_ns:21000,mad_ns:1000,min_ns:20000,max_ns:22000,checksum_verified:checksumVerified};
}
function assertRunning(ui){
  assert.ok(controlIds.every(id=>ui.get(id).disabled));assert.equal(ui.get('live-stop').disabled,false);assert.equal(ui.get('live-panel').attributes['aria-busy'],'true');assert.equal(ui.timers.size,1);
}
function assertFinished(ui,worker){
  assert.ok(controlIds.every(id=>!ui.get(id).disabled));assert.equal(ui.get('live-stop').disabled,true);assert.equal(ui.get('live-panel').attributes['aria-busy'],'false');assert.equal(worker.terminations,1);assert.equal(ui.revoked.filter(url=>url===worker.url).length,1);assert.equal(ui.timers.size,0);
}

test('a failed five-sample rotation run retries with a fresh three-sample worker and ignores stale messages',()=>{
  const ui=page();ui.choose('live-operation','rotation3d');assert.deepEqual(ui.get('live-size').options.map(option=>option.value),['3']);
  ui.choose('live-samples','5');ui.click('live-run');const first=ui.workers[0];assertRunning(ui);
  assert.deepEqual(first.messages[0].config,{operation:'rotation3d',size:3,samples:5,seed:17});assert.equal(first.messages[0].type,'benchmark');
  first.deliver({type:'progress',message:'Measured sample 2/5',check:{name:'First run partial check',passed:true}});
  assert.equal(ui.get('live-check-rows').children.length,1);
  first.deliver({type:'error',message:'browser timer resolution is too coarse for this workload',result:benchmark(5,false)});assertFinished(ui,first);
  assert.match(ui.get('live-status').textContent,/failed/i);assert.equal(ui.get('live-result-rows').children.length,0);assert.equal(ui.get('live-results').hidden,true);assert.equal(ui.get('live-export').disabled,true);

  ui.choose('live-samples','3');ui.click('live-run');const second=ui.workers[1];assertRunning(ui);
  assert.notEqual(second,first);assert.notEqual(second.url,first.url);assert.equal(ui.workers.length,2);assert.equal(second.messages.length,1);
  assert.deepEqual(second.messages[0].config,{operation:'rotation3d',size:3,samples:3,seed:17});assert.equal(first.messages[0].config.samples,5);
  assert.equal(ui.get('live-check-rows').children.length,0);assert.match(ui.get('live-check-summary').textContent,/checking/i);
  const pending=ui.snapshot();
  const deliverStale=()=>{
    first.deliver({type:'progress',message:'Stale five-sample progress',check:{name:'Stale check',passed:false}});
    first.deliver({type:'error',message:'Stale five-sample failure',checks:checks()});
    first.deliver({type:'done',kind:'benchmark',checks:checks(),result:benchmark(5)});
    first.fail('Stale worker runtime failure');
  };
  deliverStale();assert.equal(ui.snapshot(),pending);assert.equal(second.terminations,0);assert.equal(ui.timers.size,1);
  second.deliver({type:'progress',message:'Measured sample 1/3',check:{name:'Second run check',passed:true}});
  second.deliver({type:'done',kind:'benchmark',checks:checks(),result:benchmark(3)});assertFinished(ui,second);
  assert.match(ui.get('live-status').textContent,/rotation3d.*verified/i);assert.equal(ui.get('live-result-rows').children.length,1);
  assert.equal(ui.get('live-result-rows').children[0].children[5].textContent,'3 × 7');assert.equal(ui.get('live-check-rows').children.length,2);
  assert.doesNotMatch(ui.get('live-check-rows').textContent,/First run|Stale/);
  const complete=ui.snapshot();deliverStale();assert.equal(ui.snapshot(),complete);assert.equal(ui.workers.length,2);
  const exported=ui.export();assert.equal(exported.benchmarks.length,1);assert.equal(exported.benchmarks[0].operation,'rotation3d');assert.equal(exported.benchmarks[0].samples.length,3);assert.deepEqual(exported.benchmarks[0].samples,benchmark(3).samples);
  assert.equal(exported.checks.passed,2);assert.doesNotMatch(ui.get('live-export-json').textContent,/Stale|First run|browser timer resolution is too coarse for this workload/);
});

test('failed or unverified timings never enter benchmark history or exports even when correctness checks completed',()=>{
  for(const response of [{type:'error',message:'Measured checksum failed.',checks:checks(),result:benchmark(5,false)},{type:'done',kind:'benchmark',checks:checks(),result:benchmark(5,false)}]){
    const ui=page();ui.choose('live-operation','rotation3d');ui.choose('live-samples','5');ui.click('live-run');const worker=ui.workers[0];
    worker.deliver(response);assertFinished(ui,worker);assert.match(ui.get('live-status').textContent,/failed/i);
    assert.equal(ui.get('live-result-rows').children.length,0);assert.equal(ui.get('live-results').hidden,true);
    const exported=ui.export();assert.deepEqual(exported.benchmarks,[]);assert.equal(exported.checks.passed,2);
    assert.ok(!ui.get('live-export-json').textContent.includes('median_ns'));
  }
});
