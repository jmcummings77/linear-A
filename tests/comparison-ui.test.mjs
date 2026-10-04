import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';
import {test} from 'node:test';
const base=new URL('../',import.meta.url),read=p=>readFileSync(new URL(p,base),'utf8');
const report=JSON.parse(read('benchmarks/reports/wasm-profile-check/results.json'));
function harness({mobile=false}={}){
 const nodes=new Map(),downloads=[],blobs=new Map();let next=0;
 class Element{
  constructor(tag='div'){this.tagName=tag;this.children=[];this.value='';this.textContent='';this.style={};this.attributes={};this.clientWidth=800;this.width=800;this.height=180;}
  append(...children){this.children.push(...children);if(this.tagName==='select'&&this.value===''&&this.children.length)this.value=String(this.children[0].value);}
  replaceChildren(...children){this.children=[];if(this.tagName==='select')this.value='';this.append(...children);}
  get options(){return this.children;}
  setAttribute(name,value){this.attributes[name]=value;}
  focus(){this.focused=true;}
  scrollIntoView(){this.scrolled=true;}
  click(){if(this.tagName==='a')downloads.push({href:this.href,download:this.download});}
  getContext(){return new Proxy({},{get:()=>()=>{},set:()=>true});}
 }
 const selectIds=['baseline','candidate','implementation','operation','size','visibility','sort','report-theme'];
 const get=id=>{if(!nodes.has(id))nodes.set(id,new Element(selectIds.includes(id)?'select':'div'));return nodes.get(id);};
 get('snapshots').textContent=JSON.stringify([{label:'baseline',data:report},{label:'candidate',data:report}]);get('visibility').value='all';get('sort').value='identity';
 const context={document:{getElementById:get,createElement:t=>new Element(t),documentElement:{dataset:{}}},Blob,URL:{createObjectURL:b=>{const u='blob:'+ ++next;blobs.set(u,b);return u;},revokeObjectURL:u=>blobs.delete(u)},matchMedia:()=>({matches:mobile}),getComputedStyle:()=>({getPropertyValue:()=> '#888'}),devicePixelRatio:1,addEventListener(){},setTimeout(){}};
 vm.runInNewContext(read('benchmarks/comparison/model.mjs').replaceAll('export ','')+'\n'+read('benchmarks/comparison/app.mjs'),context);
 return {get,downloads,blobs};
}
test('file handler loads local JSON, invalid files preserve state, export creates numeric artifact',async()=>{
 const h=harness(),file={name:'candidate.json',size:1000,text:async()=>JSON.stringify(report)};
 await h.get('candidate-file').onchange({target:{files:[file],value:'file'}});
 assert.match(h.get('load-status').textContent,/Loaded candidate.json locally/);
 assert.match(h.get('counts').textContent,/14 comparable/);
 const previous=h.get('candidate').value;
 await h.get('candidate-file').onchange({target:{files:[{...file,text:async()=>'{broken'}],value:'file'}});
 assert.match(h.get('load-status').textContent,/File rejected/);assert.equal(h.get('candidate').value,previous);
 h.get('download').onclick();assert.equal(h.downloads.length,1);
 const data=JSON.parse(await h.blobs.get(h.downloads[0].href).text());assert.equal(data.kind,'descriptive_snapshot_comparison');assert.equal(data.rows.length,14);assert.ok(data.rows.every(r=>r.ratio===1));
});
test('local file handler enforces size limit and metadata conflicts block ratios',async()=>{
 const h=harness();await h.get('baseline-file').onchange({target:{files:[{size:21*1024*1024,name:'large.json',text:async()=>{throw Error('must not read');}}],value:'file'}});
 assert.match(h.get('load-status').textContent,/20 MiB/);
 const changed=structuredClone(report);changed.seed++;
 await h.get('candidate-file').onchange({target:{files:[{name:'changed.json',size:100,text:async()=>JSON.stringify(changed)}],value:'file'}});
 assert.match(h.get('blocked').textContent,/seeds differ/);assert.match(h.get('counts').textContent,/0 comparable/);
});
test('sample inspection follows a visible workload after filtering and clears an empty selection',()=>{
 const h=harness(),rows=h.get('rows').children;
 assert.equal(rows[0].className,'selected');
 rows[1].children[0].children[0].onclick();
 assert.equal(rows[0].className,'');assert.equal(rows[1].className,'selected');
 assert.equal(rows[1].children[0].children[0].attributes['aria-pressed'],'true');
 const operation=h.get('operation').options.at(-1).value;h.get('operation').value=operation;h.get('operation').onchange();
 assert.ok(h.get('sample-title').textContent.includes(operation));
 assert.equal(h.get('rows').children.filter(row=>row.className==='selected').length,1);
 h.get('operation').value='no recorded operation';h.get('operation').onchange();
 assert.equal(h.get('rows').children.length,0);assert.match(h.get('sample-title').textContent,/No workloads match/);assert.equal(h.get('sample-values').textContent,'');
});
test('mobile workload activation moves to the inspector and provides a return link',()=>{
 const h=harness({mobile:true});assert.equal(h.get('sample-panel').scrolled,undefined);
 const button=h.get('rows').children[1].children[0].children[0];button.onclick();
 assert.equal(h.get('sample-panel').scrolled,true);assert.equal(h.get('sample-heading').focused,true);
 assert.equal(h.get('sample-back').href,'#'+button.id);
});
test('snapshot identities follow loaded files and swap consistently in results and inspector',async()=>{
 const h=harness(),changed=structuredClone(report);changed.created_at='2026-10-04T01:02:03Z';changed.revision='abcdef0123456789abcdef0123456789abcdef01';changed.dirty=true;
 const label='<img src=x> candidate.json';
 await h.get('candidate-file').onchange({target:{files:[{name:label,size:100,text:async()=>JSON.stringify(changed)}],value:'file'}});
 for(const prefix of ['', 'sample-']){
  assert.equal(h.get(prefix+'candidate-label').textContent,label);assert.equal(h.get(prefix+'candidate-label').children.length,0);
  assert.equal(h.get(prefix+'candidate-reference').textContent,'2026-10-04 · abcdef01 · dirty');
 }
 h.get('swap').onclick();
 for(const prefix of ['', 'sample-']){
  assert.equal(h.get(prefix+'baseline-label').textContent,label);
  assert.equal(h.get(prefix+'baseline-reference').textContent,'2026-10-04 · abcdef01 · dirty');
  assert.equal(h.get(prefix+'candidate-label').textContent,'baseline');
 }
});
