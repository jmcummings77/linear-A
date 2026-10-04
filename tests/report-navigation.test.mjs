import assert from 'node:assert/strict';
import test from 'node:test';
import vm from 'node:vm';
import {readFileSync} from 'node:fs';
const source=readFileSync(new URL('../benchmarks/report-navigation.js',import.meta.url),'utf8');
const css=readFileSync(new URL('../benchmarks/report-navigation.css',import.meta.url),'utf8');

// Minimal DOM fixture implements native node identity, parentage, attributes and anchor traversal.
class Element {
 constructor(tag='div',id='',text=''){this.tagName=tag.toUpperCase();this.attributes=new Map();this.children=[];this.parentElement=null;this.ownText=text;if(id)this.id=id;}
 get id(){return this.getAttribute('id')||'';} set id(v){this.setAttribute('id',v);}
 get className(){return this.getAttribute('class')||'';} set className(v){this.setAttribute('class',v);}
 get classList(){const node=this;return {contains:c=>node.className.split(/\s+/).includes(c),add:c=>{if(!node.classList.contains(c))node.className=[node.className,c].filter(Boolean).join(' ');},toggle(c,on){if(on)node.classList.add(c);else node.className=node.className.split(/\s+/).filter(v=>v!==c).join(' ');}};}
 setAttribute(k,v){this.attributes.set(k,String(v));} getAttribute(k){return this.attributes.get(k)??null;} removeAttribute(k){this.attributes.delete(k);}
 get hidden(){return this.attributes.has('hidden');}set hidden(v){v?this.setAttribute('hidden',''):this.removeAttribute('hidden');}
 get textContent(){return this.ownText+this.children.map(node=>node.textContent).join('');}set textContent(v){this.ownText=v;this.replaceChildren();}
 detach(){if(this.parentElement)this.parentElement.children=this.parentElement.children.filter(node=>node!==this);this.parentElement=null;}
 append(...nodes){for(const node of nodes){node.detach();this.children.push(node);node.parentElement=this;}}
 replaceChildren(...nodes){for(const child of this.children)child.parentElement=null;this.children=[];this.append(...nodes);}
 replaceWith(...nodes){const parent=this.parentElement,index=parent.children.indexOf(this);this.detach();for(const node of nodes)node.detach();parent.children.splice(index,0,...nodes);for(const node of nodes)node.parentElement=parent;}
 matches(selector){
  if(selector.includes(','))return selector.split(',').some(s=>this.matches(s.trim()));
  if(selector.startsWith('.'))return this.classList.contains(selector.slice(1));
  if(selector.startsWith('#'))return this.id===selector.slice(1);
  if(selector.startsWith('['))return this.attributes.has(selector.slice(1,-1));
  if(selector==='a[href]')return this.tagName==='A'&&this.attributes.has('href');
  return this.tagName===selector.toUpperCase();
 }
 querySelectorAll(selector){
  const split=selector.lastIndexOf(' ');
  if(split>0&&!selector.includes(',')){const ancestor=selector.slice(0,split),last=selector.slice(split+1);return this.querySelectorAll(last).filter(n=>n.parentElement?.closest(ancestor));}
  return this.children.flatMap(node=>[...(node.matches(selector)?[node]:[]),...node.querySelectorAll(selector)]);
 }
 querySelector(s){return this.querySelectorAll(s)[0]||null;}
 closest(s){for(let n=this;n;n=n.parentElement)if(n.matches(s))return n;return null;}
 contains(other){for(let n=other;n;n=n.parentElement)if(n===this)return true;return false;}
 scrollIntoView(){this.scrollCount=(this.scrollCount||0)+1;}
}
function fixture({family='core',hash='',hidden=[],ready='complete',directory=false}={}){
 const root=new Element('body'),main=new Element('main'),header=new Element('header','overview','Identity'),topbar=new Element('div');topbar.className='report-topbar';
 const link=(href,text)=>{const n=new Element('a','',text);n.setAttribute('href',href);return n;};
 const back=link('../','All reports');topbar.append(back);root.append(main);main.append(topbar,header);
 const nav=new Element('nav');nav.className='report-nav';main.append(nav);
 const ids=family==='core'?['timings-panel','profiles-panel','correctness-panel','method-panel','explore','live-panel','accuracy-panel','geometry-panel']:family==='sparse'?['runtime','reuse-section','storage-section','method-section','explore']:family==='plain'?['snapshots-panel','results','methodology']:['results','explore','method-section'];
 if(directory)main.append(new Element('div','report-library'));
 const get=id=>root.querySelector('#'+id);
 for(const id of ids){const section=new Element('section',id);section.append(new Element('h2','',id));section.hidden=hidden.includes(id);main.append(section);nav.append(link('#'+id,id));}
 if(family==='core'){const stats=new Element('div','stats');main.append(stats);const context=new Element('div','saved-context');main.append(context);const note=new Element('p','','Availability');note.className='availability-note';main.append(note);const local=new Element('nav');local.append(link('#geometry-panel','Geometry'));get('explore').append(local);}
 if(family==='sparse'){header.id='';const summary=new Element('div','overview');main.append(summary);const label=new Element('div');label.className='explore-label';main.append(label);}
 const download=link('results.json','Download JSON');download.setAttribute('download','');nav.append(download);
 const related=new Element('nav');related.className='related-reports';const sibling=link('../ic0/','IC(0)');related.append(sibling);main.append(related);
 const filter=new Element('select','filter');filter.value='rust';get(ids[0]).append(filter);
 const details=new Element('details');const nested=new Element('pre','provenance-detail','Source');details.append(nested);get(family==='core'?'method-panel':family==='plain'?'methodology':'method-section').append(details);
 const location=new URL('https://example.test/latest/'+hash),listeners=new Map(),documentListeners=new Map(),frames=new Map(),observers=[];let nextFrame=0,resizes=0;
 const add=(map,k,fn)=>map.set(k,[...(map.get(k)||[]),fn]);
 const emit=(map,k,event={})=>{for(const fn of map.get(k)||[])fn(event);};
 const window={location,addEventListener:(k,f)=>add(listeners,k,f),requestAnimationFrame:f=>{frames.set(++nextFrame,f);return nextFrame;},cancelAnimationFrame:id=>frames.delete(id),dispatchEvent:event=>{if(event.type==='resize')resizes++;emit(listeners,event.type,event);}};
 const document={readyState:ready,querySelector:s=>root.querySelector(s),getElementById:get,createElement:t=>new Element(t),addEventListener:(k,f)=>add(documentListeners,k,f)};
 class MutationObserver{constructor(callback){this.callback=callback;observers.push(this);}observe(){}}
 const context={document,window,URL,Event,MutationObserver};
 const boot=()=>vm.runInNewContext(source,context);boot();
 const flush=()=>{const queued=[...frames.values()];frames.clear();queued.forEach(f=>f());};
 const navigate=(hash,event='hashchange')=>{location.hash=hash;emit(listeners,event);flush();};
 const click=(anchor,options={})=>{emit(documentListeners,'click',{target:anchor,button:0,...options});};
 return {main,header,nav,back,filter,download,sibling,related,get,details,nested,window,document,observers,flush,navigate,click,boot,link,
  ready:()=>{document.readyState='interactive';emit(documentListeners,'DOMContentLoaded');},
  active:()=>main.getAttribute('data-report-view'),modeNav:()=>main.querySelector('.report-view-nav'),localNav:()=>main.querySelector('.report-view-sections'),
  get resizes(){return resizes;},isInactive:id=>get(id).classList.contains('report-mode-inactive')};
}

test('core defaults to results and retains the header, controls and native hidden states',()=>{
 const h=fixture({hidden:['accuracy-panel']});assert.equal(h.active(),'results');assert.ok(!h.isInactive('timings-panel'));assert.ok(h.isInactive('geometry-panel'));assert.ok(h.isInactive('method-panel'));assert.ok(!h.header.classList.contains('report-mode-inactive'));assert.equal(h.get('accuracy-panel').hidden,true);
 assert.equal(h.get('filter'),h.filter);assert.equal(h.filter.value,'rust');assert.equal(h.modeNav().getAttribute('aria-label'),'Report views');assert.equal(h.modeNav().getAttribute('role'),null);assert.equal(h.modeNav().children.filter(n=>n.getAttribute('aria-current')==='page').length,1);
 assert.equal(h.back.textContent,'Report library');assert.equal(h.back.getAttribute('href'),'../');assert.equal(h.back.getAttribute('data-report-icon'),'back');
});
test('old demo and nested method deep links reveal their view and preserve the URL',()=>{
 const h=fixture({hash:'#geometry-panel'});h.flush();assert.equal(h.active(),'demo');assert.ok(h.get('geometry-panel').scrollCount);assert.equal(h.window.location.hash,'#geometry-panel');
 h.navigate('#provenance-detail');assert.equal(h.active(),'evidence');assert.equal(h.details.open,true);assert.ok(h.nested.scrollCount);assert.equal(h.window.location.hash,'#provenance-detail');
 h.navigate('#explore');assert.equal(h.active(),'demo');
});
test('native link click reveals before scrolling and hash/popstate restore back and forward views',()=>{
 const h=fixture();const demo=h.modeNav().children[1];h.click(demo);assert.equal(h.active(),'demo');assert.equal(h.window.location.hash,'');
 h.navigate(demo.getAttribute('href'));h.navigate('#method-panel');assert.equal(h.active(),'evidence');h.navigate('#live-panel','popstate');assert.equal(h.active(),'demo');h.navigate('','popstate');assert.equal(h.active(),'results');
 assert.equal(h.get('filter'),h.filter);assert.equal(h.filter.value,'rust');
});
test('modified clicks do not change the current view',()=>{
 const h=fixture();h.click(h.modeNav().children[1],{ctrlKey:true});assert.equal(h.active(),'results');h.click(h.modeNav().children[1],{button:1});assert.equal(h.active(),'results');
});
test('verification-only captures omit empty results/demo views without revealing unavailable sections',()=>{
 const h=fixture({hidden:['timings-panel','profiles-panel','explore','live-panel','accuracy-panel','geometry-panel']});assert.equal(h.active(),'evidence');assert.deepEqual(h.modeNav().children.map(n=>n.textContent),['Method & data']);assert.equal(h.get('timings-panel').hidden,true);assert.equal(h.get('geometry-panel').hidden,true);assert.ok(!h.isInactive('stats'));
 h.navigate('#geometry-panel');assert.equal(h.active(),'evidence');assert.equal(h.get('geometry-panel').scrollCount,undefined);
});
test('local anchors appear only for multiple available sections and track nested targets',()=>{
 const h=fixture({hidden:['profiles-panel','accuracy-panel']});assert.equal(h.localNav().hidden,true);h.navigate('#geometry-panel');assert.equal(h.localNav().hidden,false);const anchors=h.localNav().querySelectorAll('a[href]');assert.deepEqual(anchors.map(n=>n.getAttribute('href')),['#live-panel','#geometry-panel']);assert.equal(anchors[1].getAttribute('aria-current'),'location');assert.ok(h.get('explore').querySelector('nav').classList.contains('report-view-redundant-nav'));
});
test('sparse dynamic availability updates local links without resetting filters or other hidden states',()=>{
 const h=fixture({family:'sparse',hidden:['reuse-section']});assert.equal(h.active(),'results');assert.equal(h.localNav().querySelectorAll('a[href]').length,2);h.get('reuse-section').hidden=false;h.observers[0].callback();assert.equal(h.localNav().querySelectorAll('a[href]').length,3);assert.equal(h.filter.value,'rust');h.navigate('#explore');assert.ok(h.isInactive('overview'));assert.equal(h.localNav().hidden,true);
});
test('ordering, Cholesky, AMD, IC0 and multigrid section pattern has three nonempty views',()=>{
 const h=fixture({family:'special'});assert.deepEqual(h.modeNav().children.map(n=>n.textContent),['Results','Interactive demo','Method & data']);assert.equal(h.localNav().hidden,true);h.navigate('#method-section');assert.equal(h.active(),'evidence');assert.ok(!h.isInactive('method-section'));assert.ok(h.isInactive('results'));
});
test('downloads and related destinations are outside the mode bar and preserve their links',()=>{
 const h=fixture({family:'special'});assert.equal(h.download.parentElement.className,'report-view-utilities');assert.equal(h.sibling.getAttribute('href'),'../ic0/');assert.equal(h.related.parentElement.tagName,'DETAILS');assert.equal(h.related.parentElement.querySelector('summary').textContent,'Related reports');assert.equal(h.modeNav().children.length,3);assert.equal(h.nav.parentElement,null);
});
test('revealing a new view requests chart resizing; profile shortcut still works',()=>{
 const h=fixture();h.flush();const before=h.resizes;h.navigate('#geometry-panel');assert.ok(h.resizes>before);h.window.dispatchEvent(new Event('matrix-profile-focus'));assert.equal(h.active(),'results');assert.equal(h.window.location.hash,'#profiles-panel');h.flush();assert.ok(h.resizes>before+1);
});
test('standalone comparisons keep all sections and gain only active anchor state',()=>{
 const h=fixture({family:'plain',hash:'#results'});assert.equal(h.active(),null);assert.equal(h.nav.parentElement,h.main);assert.equal(h.nav.children.find(n=>n.getAttribute('href')==='#results').getAttribute('aria-current'),'location');h.navigate('#methodology');assert.equal(h.nav.children.find(n=>n.getAttribute('href')==='#methodology').getAttribute('aria-current'),'location');assert.ok(!h.main.classList.contains('report-modes-enabled'));
});
test('boot waits for existing scripts, is idempotent, and CSS leaves unenhanced pages available',()=>{
 const h=fixture({ready:'loading'});assert.equal(h.active(),null);assert.ok(!h.get('timings-panel').classList.contains('report-mode-inactive'));h.ready();assert.equal(h.active(),'results');h.boot();assert.equal(h.main.querySelectorAll('.report-view-nav').length,1);assert.match(css,/main\.report-modes-enabled > \.report-mode-inactive/);assert.match(css,/@media print/);assert.match(css,/min-height:44px/);
});

test('the question-first directory never receives modes or a report back label',()=>{
 const h=fixture({directory:true});assert.equal(h.active(),null);assert.equal(h.nav.parentElement,h.main);assert.equal(h.back.textContent,'All reports');assert.ok(!h.main.classList.contains('report-modes-enabled'));
});
test('standalone reports retain their original route with the common library label',()=>{
 const h=fixture({family:'plain'});assert.equal(h.back.textContent,'Report library');assert.equal(h.back.getAttribute('href'),'../');assert.equal(h.back.getAttribute('data-report-icon'),'back');
});
test('view links have visible names and consistent decorative action icons',()=>{
 const h=fixture();assert.deepEqual(h.modeNav().children.map(n=>n.getAttribute('data-report-icon')),['chart','play','info']);assert.ok(h.modeNav().children.every(n=>n.textContent.length));
});
