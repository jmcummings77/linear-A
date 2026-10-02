import assert from 'node:assert/strict';
import test from 'node:test';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';

const source=readFileSync(new URL('../benchmarks/timing-report.js',import.meta.url),'utf8');

/** Small DOM model with native select selection and ordinary event handlers. */
class Element {
  constructor(tag){this.tagName=tag.toUpperCase();this.children=[];this.parentElement=null;this.dataset={};this.style={};this.className='';this.attributes={};this.listeners=new Map();this._text='';this._value='';this.selectedIndex=-1;this.disabled=false;this.checked=false;this.open=false;}
  set textContent(value){this.children.forEach(x=>{x.parentElement=null;});this.children=[];this._text=String(value);}
  get textContent(){return this._text+this.children.map(x=>x.textContent).join('');}
  set value(value){const text=String(value);if(this.tagName==='SELECT')this.selectedIndex=this.options.findIndex(x=>x.value===text);else this._value=text;}
  get value(){return this.tagName==='SELECT'?(this.options[this.selectedIndex]?.value??''):this._value;}
  get options(){return this.children.filter(x=>x.tagName==='OPTION');}
  append(...nodes){for(let node of nodes){if(typeof node==='string'){const text=new Element('#text');text.textContent=node;node=text;}node.parentElement=this;this.children.push(node);if(this.tagName==='SELECT'&&this.selectedIndex<0)this.selectedIndex=0;}}
  appendChild(node){this.append(node);return node;}
  replaceChildren(...nodes){this.children.forEach(x=>{x.parentElement=null;});this.children=[];this._text='';if(this.tagName==='SELECT')this.selectedIndex=-1;this.append(...nodes);}
  setAttribute(name,value){this.attributes[name]=String(value);if(name==='class')this.className=String(value);}
  getAttribute(name){return this.attributes[name]??null;}
  contains(node){return this===node||this.children.some(x=>x.contains(node));}
  addEventListener(type,handler){if(!this.listeners.has(type))this.listeners.set(type,[]);this.listeners.get(type).push(handler);}
  dispatch(type,extra={}){const event={type,target:this,preventDefault(){},stopPropagation(){},...extra};this['on'+type]?.(event);for(const listener of this.listeners.get(type)||[])listener(event);return event;}
  click(){this.dispatch('click');}
  focus(){this.focused=true;}
  querySelectorAll(selector){const result=[];const matches=node=>selector==='input'||selector==='input[type="checkbox"]'?node.tagName==='INPUT'&&(selector==='input'||node.type==='checkbox'):selector==='[data-language]'?node.dataset.language!==undefined:selector.startsWith('.')?node.className.split(' ').includes(selector.slice(1)):node.tagName===selector.toUpperCase();for(const child of this.children){if(matches(child))result.push(child);result.push(...child.querySelectorAll(selector));}return result;}
  querySelector(selector){return this.querySelectorAll(selector)[0]??null;}
}

const timing=(implementation,operation,size,median)=>({implementation,operation,size,status:'passed',median_ns:median,mad_ns:median/10,min_ns:median*.9,max_ns:median*1.1,iterations:10,samples:[{ns_per_op:median*.9},{ns_per_op:median*1.1}]});
function fixture(){return {implementations:[{id:'csharp',name:'C#',status:'passed'},{id:'rust',name:'Rust',status:'passed'},{id:'python',name:'Python',status:'passed'},{id:'java',name:'Java',status:'failed'}],results:[timing('csharp','multiply',8,10),timing('rust','multiply',8,20),timing('python','multiply',8,100),timing('csharp','multiply',16,30),timing('rust','multiply',16,60),timing('python','multiply',16,300),timing('csharp','transpose',8,2),timing('python','transpose',8,10),{...timing('java','multiply',8,1000),status:'failed'}]};}
function report(data=fixture()){
  const original=JSON.stringify(data),nodes=new Map(),listeners=new Map();
  for(const [id,tag]of Object.entries({operation:'select',size:'select',baseline:'select',chart:'div',samples:'table','chart-note':'p','timing-languages':'details','timing-language-summary':'summary','timing-language-options':'fieldset','timing-language-all':'button','timing-language-none':'button',profile:'select','profile-note':'p'}))nodes.set(id,new Element(tag));
  const get=id=>{assert.ok(nodes.has(id),`unknown DOM id ${id}`);return nodes.get(id);};
  get('timing-languages').append(get('timing-language-summary'),get('timing-language-options'),get('timing-language-all'),get('timing-language-none'));
  get('profile-note').textContent='Untouched profile selection';
  const document={getElementById:get,createElement:tag=>new Element(tag),createTextNode:text=>{const node=new Element('#text');node.textContent=text;return node;},addEventListener(type,handler){if(!listeners.has(type))listeners.set(type,[]);listeners.get(type).push(handler);}};
  const el=(tag,text,cls)=>{const node=document.createElement(tag);if(text!==undefined)node.textContent=text;if(cls)node.className=cls;return node;};
  const addOption=(select,value,text)=>{const option=el('option',text);option.value=value;select.append(option);};
  const context={data,passed:data.implementations.filter(x=>x.status==='passed'),valid:data.results.filter(x=>x.status==='passed'),names:Object.fromEntries(data.implementations.map(x=>[x.id,x.name])),$:get,el,addOption,time:n=>`${n.toFixed(1)} ns`,speed:n=>n.toFixed(2),document};
  vm.runInNewContext(source,context);
  const checkboxes=()=>get('timing-language-options').querySelectorAll('input');
  return {get,data,original,checkboxes,
    choose(id,value){get(id).value=value;get(id).dispatch('change');},
    language(id,checked){const input=checkboxes().find(x=>x.dataset.language===id);assert.ok(input,`missing language ${id}`);input.checked=checked;input.dispatch('change');},
    click(id){get(id).click();},
    event(type,event){for(const callback of listeners.get(type)||[])callback({type,preventDefault(){},...event});},
    rows(){return get('chart').children.filter(x=>x.className.split(' ').includes('bar-row'));},
    labels(){return this.rows().map(x=>x.children[0].textContent);},
    widths(){return this.rows().map(row=>row.querySelector('.bar').style.width);},
    baseline(){return get('baseline').options.map(x=>x.value);},
    sampleLabels(){return get('samples').children.slice(1).map(row=>row.children[0].textContent);},
  };
}

test('timing languages start checked, show all successful rows, and leave recorded data intact',()=>{
  const page=report();
  assert.deepEqual(page.checkboxes().map(x=>x.dataset.language).sort(),['csharp','python','rust']);
  assert.ok(page.checkboxes().every(x=>x.type==='checkbox'&&x.checked));
  assert.deepEqual(page.labels(),['C#','Rust','Python']);
  assert.deepEqual(page.widths(),['10%','20%','100%']);
  assert.deepEqual(page.sampleLabels(),page.labels());
  assert.equal(page.get('baseline').value,'csharp');
  assert.equal(page.get('baseline').disabled,false);
  assert.match(page.get('timing-language-summary').textContent,/3/);
  assert.equal(JSON.stringify(page.data),page.original);
});

test('unchecking a slow language hides its chart, samples and baseline and rescales visible bars',()=>{
  const page=report();page.choose('baseline','python');page.language('python',false);
  assert.deepEqual(page.labels(),['C#','Rust']);assert.deepEqual(page.sampleLabels(),page.labels());
  assert.deepEqual(page.widths(),['50%','100%']);
  assert.deepEqual(page.baseline(),['csharp','rust']);assert.equal(page.get('baseline').value,'csharp');
  assert.match(page.get('timing-language-summary').textContent,/2/);
  assert.equal(page.get('profile-note').textContent,'Untouched profile selection');
  assert.equal(JSON.stringify(page.data),page.original);
});

test('language selections persist across sizes and operations and visible baselines survive',()=>{
  const page=report();page.language('python',false);page.choose('baseline','rust');
  page.choose('size','16');assert.deepEqual(page.labels(),['C#','Rust']);assert.equal(page.get('baseline').value,'rust');
  page.choose('operation','transpose');assert.deepEqual(page.labels(),['C#']);assert.equal(page.get('baseline').value,'csharp');
  assert.equal(page.checkboxes().find(x=>x.dataset.language==='python').checked,false);
  page.choose('operation','multiply');assert.deepEqual(page.labels(),['C#','Rust']);assert.deepEqual(page.widths(),['50%','100%']);
  page.language('python',true);assert.deepEqual(page.labels(),['C#','Rust','Python']);assert.deepEqual(page.widths(),['10%','20%','100%']);
});

test('none empties only timing results, disables its baseline and clears stale notes; all restores results',()=>{
  const page=report();assert.ok(page.get('chart-note').textContent.length>0);page.click('timing-language-none');
  assert.ok(page.checkboxes().every(x=>!x.checked));assert.deepEqual(page.rows(),[]);assert.equal(page.get('samples').children.length,0);
  assert.equal(page.get('baseline').disabled,true);assert.deepEqual(page.baseline(),[]);assert.equal(page.get('chart-note').textContent,'');
  assert.match(page.get('chart').textContent,/select|language/i);assert.match(page.get('timing-language-summary').textContent,/0|none/i);
  page.choose('size','16');assert.deepEqual(page.rows(),[]);assert.equal(page.get('baseline').disabled,true);
  page.click('timing-language-all');assert.ok(page.checkboxes().every(x=>x.checked));assert.equal(page.get('baseline').disabled,false);
  assert.deepEqual(page.labels(),['C#','Rust','Python']);assert.deepEqual(page.widths(),['10%','20%','100%']);
  assert.equal(page.get('profile-note').textContent,'Untouched profile selection');assert.equal(JSON.stringify(page.data),page.original);
});

test('one language checkbox groups only known determinant suffixes across algorithm variants',()=>{
  const ids=['go-lu','go-cholesky','go-cofactor','wasm','go-simd'];
  const names=['Go / LU','Go / Cholesky','Go / cofactor','WebAssembly','Go SIMD'];
  const page=report({implementations:ids.map((id,i)=>({id,name:names[i],status:'passed'})),results:ids.map((id,i)=>timing(id,'determinant',8,10*(i+1)))});
  assert.deepEqual(page.checkboxes().map(x=>x.dataset.language).sort(),['go','go-simd','wasm']);
  assert.equal(page.checkboxes().find(x=>x.dataset.language==='go').parentElement.textContent,'Go');
  page.choose('baseline','go-cholesky');page.language('go',false);
  assert.deepEqual(page.labels(),['WebAssembly','Go SIMD']);assert.deepEqual(page.sampleLabels(),page.labels());
  assert.deepEqual(page.baseline(),['wasm','go-simd']);assert.deepEqual(page.widths(),['80%','100%']);
  page.language('go',true);assert.deepEqual(page.labels(),names);assert.equal(page.checkboxes().length,3);
});

test('zero-duration rows have finite zero-width bars and never display a NaN ratio',()=>{
  const data=fixture();data.results=[timing('csharp','multiply',8,0),timing('rust','multiply',8,0)];
  const page=report(data);assert.deepEqual(page.widths(),['0%','0%']);assert.doesNotMatch(page.get('chart').textContent,/NaN|Infinity/);
  page.language('csharp',false);assert.deepEqual(page.widths(),['0%']);assert.equal(page.get('baseline').value,'rust');
});

test('language dropdown stays open inside and closes on outside clicks or Escape',()=>{
  const page=report(),details=page.get('timing-languages');details.open=true;
  page.event('click',{target:page.checkboxes()[0]});assert.equal(details.open,true);
  page.event('click',{target:page.get('chart')});assert.equal(details.open,false);
  details.open=true;page.event('keydown',{key:'ArrowDown',target:page.checkboxes()[0]});assert.equal(details.open,true);
  page.event('keydown',{key:'Escape',target:page.checkboxes()[0]});assert.equal(details.open,false);
});

test('missing successful timing data retains an explicit empty chart and usable language controls',()=>{
  const data=fixture();data.results=data.results.map(x=>({...x,status:'failed'}));const page=report(data);
  assert.deepEqual(page.rows(),[]);assert.equal(page.get('baseline').disabled,true);assert.equal(page.get('chart-note').textContent,'');
  assert.match(page.get('chart').textContent,/no|empty|unavailable/i);page.click('timing-language-all');page.click('timing-language-none');
  assert.deepEqual(page.rows(),[]);assert.equal(JSON.stringify(page.data),page.original);
});
