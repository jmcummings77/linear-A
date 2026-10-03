import assert from 'node:assert/strict';
import test from 'node:test';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';

const template=readFileSync(new URL('../benchmarks/report.html',import.meta.url),'utf8');
const source=template.slice(template.indexOf('let profileFrameRows='),template.indexOf('function renderFlames(){'));
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
