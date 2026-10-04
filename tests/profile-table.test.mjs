import assert from 'node:assert/strict';
import test from 'node:test';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';

const template=readFileSync(new URL('../benchmarks/report.html',import.meta.url),'utf8');
const source=template.slice(template.indexOf('let profileFrameRows='),template.indexOf('function renderFlames('));
class Element {
  constructor(){this.children=[];this._text='';this.value='';this.attributes={};this.hidden=false;this.disabled=false;}
  set textContent(value){this._text=String(value);this.children=[];}
  get textContent(){return this._text+this.children.map(child=>child.textContent).join('');}
  append(...children){this.children.push(...children);}
  replaceChildren(...children){this.children=children;this._text='';}
  setAttribute(name,value){this.attributes[name]=value;}
  focus(){this.focused=true;}
}
function harness(){
  const nodes=new Map(['profile-search','profile-frames','profile-frame-count','profile-more','reset'].map(id=>[id,new Element()]));
  const get=id=>{assert.ok(nodes.has(id),'Unknown ID '+id);return nodes.get(id);};
  const context={zoom:[],$:get,el(tag,text,cls){const node=new Element();node.tag=tag;node.className=cls;if(text!==undefined)node.textContent=text;return node;},renderFlames(){context.renders=(context.renders||0)+1;}};
  vm.createContext(context);vm.runInContext(source,context);
  return {get,context,update(stacks,total){context.updateProfileFrames(stacks,total);},rows(){return get('profile-frames').children;},search(value){get('profile-search').value=value;get('profile-search').oninput();}};
}

test('frame table aggregates inclusive weights by exact call path with explicit share',()=>{
  const page=harness();page.update([{stack:['root','multiply'],weight:3},{stack:['root','multiply'],weight:2},{stack:['root','other'],weight:5}],10);
  assert.equal(page.rows().length,3);
  const root=page.rows()[0];assert.equal(root.children[0].textContent,'root');assert.equal(root.children[1].textContent,'10.000 ms');assert.equal(root.children[2].textContent,'100.0%');
  page.search('multiply');assert.equal(page.rows().length,1);assert.match(page.rows()[0].textContent,/5.000 ms50.0%/);
  assert.equal(page.rows()[0].children[3].children[0].attributes['aria-label'],'Zoom to multiply');
});

test('ordinary button activation zooms to the frame and returns focus to reset',()=>{
  const page=harness();page.update([{stack:['root','multiply'],weight:5}],5);page.search('multiply');
  page.rows()[0].children[3].children[0].onclick();
  assert.deepEqual(Array.from(page.context.zoom),['root','multiply']);assert.equal(page.context.renders,1);assert.equal(page.get('reset').focused,true);
  page.update([{stack:['root','multiply'],weight:5}],5);
  assert.equal(page.rows()[0].children[3].children[0].disabled,true);assert.match(page.rows()[0].textContent,/Current root/);
});

test('frame search includes parent source paths, handles no matches, and paginates without dropping frames',()=>{
  const page=harness();page.update(Array.from({length:120},(_,i)=>({stack:['root (source.c:1)','child '+i],weight:1})),120);
  assert.equal(page.rows().length,100);assert.equal(page.get('profile-more').hidden,false);
  page.get('profile-more').onclick();assert.equal(page.rows().length,121);assert.equal(page.get('profile-more').hidden,true);
  page.search('source.c');assert.equal(page.rows().length,100);assert.match(page.get('profile-frame-count').textContent,/100 of 121/);
  page.search('missing-function');assert.equal(page.rows().length,0);assert.match(page.get('profile-frame-count').textContent,/0 of 0/);
});

function graphHarness(){
 const ids=['flames','profile-viewport','profile-empty','profile-note','profile','view','reset','profile-search','profile-frames','profile-frame-count','profile-more','tooltip'],nodes=new Map(ids.map(id=>[id,new Element()])),events=new Map(),get=id=>nodes.get(id);
 for(const node of nodes.values())node.style={};
 const canvas=get('flames'),viewport=get('profile-viewport');
 viewport.clientWidth=800;viewport.clientHeight=420;viewport.scrollTop=0;
 Object.defineProperty(viewport,'scrollHeight',{get:()=>parseFloat(canvas.style.height)||0});
 viewport.scrollIntoView=options=>{viewport.lastScroll=options;};
 canvas.parentElement=viewport;canvas.getContext=()=>new Proxy({},{get:()=>()=>{},set:()=>true});
 canvas.getBoundingClientRect=()=>({left:0,top:-viewport.scrollTop});
 canvas.scrollIntoView=()=>assert.fail('Scroll the bounded viewport, not the full-height canvas');
 get('profile').value='c';get('view').value='flame';get('view').options=[{},{}];
 const profiles=[{implementation:'c',chronological:true,status:'available',stacks:[Array.from({length:40},(_,i)=>i?'frame '+i:'root')],weights:[10]},{implementation:'wasm',chronological:true,status:'available',stacks:[['root','wm_multiply',...Array.from({length:30},(_,i)=>'wasm '+i)]],weights:[5]}];
 const context={profiles,data:{profiles},document:{documentElement:{}},$:get,devicePixelRatio:1,innerWidth:1000,innerHeight:800,el(tag,text,cls){const node=new Element();node.tag=tag;node.className=cls;if(text!==undefined)node.textContent=text;return node;},getComputedStyle:()=>({getPropertyValue:name=>name==='--flame-lightness'?'70':name==='--flame-saturation'?'60':'#777'}),addEventListener:(name,handler)=>events.set(name,handler)};
 vm.createContext(context);vm.runInContext(template.slice(template.indexOf('let zoom=[],rects=[];'),template.indexOf('const sectionLinks=')),context);
 return {get,viewport,canvas,context,event:name=>events.get(name)(),rects:()=>vm.runInContext('rects',context),search(value){get('profile-search').value=value;get('profile-search').oninput();}};
}

test('deep flame roots start visible and remain visible after reset, profile changes and frame zoom',()=>{
 const page=graphHarness(),viewport=page.viewport;
 const rootVisible=()=>{const root=page.rects()[0];assert.ok(root.y>=viewport.scrollTop&&root.y+root.h<=viewport.scrollTop+viewport.clientHeight,'Flame root is outside the viewport');};
 assert.ok(viewport.scrollHeight>viewport.clientHeight);rootVisible();
 viewport.scrollTop=0;page.get('reset').onclick();rootVisible();
 page.search('frame 4');page.get('profile-frames').children[0].children[3].children[0].onclick();rootVisible();
 page.get('profile').value='wasm';page.get('profile').onchange();rootVisible();
 const root=page.rects()[0];page.canvas.onclick({clientX:root.x+10,clientY:root.y-viewport.scrollTop+5});rootVisible();assert.equal(viewport.lastScroll.block,'nearest');
 page.event('matrix-profile-focus');rootVisible();assert.equal(viewport.lastScroll.block,'center');
});

test('timeline starts at top while resize and appearance changes preserve deliberate scrolling',()=>{
 const page=graphHarness(),viewport=page.viewport;
 viewport.scrollTop=125;page.event('resize');assert.equal(viewport.scrollTop,125);
 page.event('matrix-theme-change');assert.equal(viewport.scrollTop,125);
 page.get('view').value='timeline';page.get('view').onchange();assert.equal(viewport.scrollTop,0);
 viewport.scrollTop=85;page.event('matrix-theme-change');assert.equal(viewport.scrollTop,85);
 page.get('reset').onclick();assert.equal(viewport.scrollTop,0);
 page.get('view').value='flame';page.get('view').onchange();assert.equal(viewport.scrollTop,viewport.scrollHeight-viewport.clientHeight);
 viewport.clientHeight=360;page.event('resize');assert.ok(viewport.scrollTop>=0&&viewport.scrollTop<=viewport.scrollHeight-viewport.clientHeight);
});
