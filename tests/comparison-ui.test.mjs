import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';
import {test} from 'node:test';
const base=new URL('../',import.meta.url),read=p=>readFileSync(new URL(p,base),'utf8');
const report=JSON.parse(read('benchmarks/reports/wasm-profile-check/results.json'));
function harness(){
 const nodes=new Map(),downloads=[],blobs=new Map();let next=0;
 class Element{
  constructor(tag='div'){this.tagName=tag;this.children=[];this.value='';this.textContent='';this.style={};this.clientWidth=800;this.width=800;this.height=180;}
  append(...children){this.children.push(...children);if(this.tagName==='select'&&this.value===''&&this.children.length)this.value=String(this.children[0].value);}
  replaceChildren(...children){this.children=[];if(this.tagName==='select')this.value='';this.append(...children);}
  get options(){return this.children;}
  click(){if(this.tagName==='a')downloads.push({href:this.href,download:this.download});}
  getContext(){return new Proxy({},{get:()=>()=>{},set:()=>true});}
 }
 const selectIds=['baseline','candidate','implementation','operation','size','visibility','sort','theme'];
 const get=id=>{if(!nodes.has(id))nodes.set(id,new Element(selectIds.includes(id)?'select':'div'));return nodes.get(id);};
 get('snapshots').textContent=JSON.stringify([{label:'baseline',data:report},{label:'candidate',data:report}]);get('visibility').value='all';get('sort').value='identity';get('theme').value='system';
 const context={document:{getElementById:get,createElement:t=>new Element(t),documentElement:{dataset:{}}},Blob,URL:{createObjectURL:b=>{const u='blob:'+ ++next;blobs.set(u,b);return u;},revokeObjectURL:u=>blobs.delete(u)},localStorage:{getItem:()=>null,setItem(){}},matchMedia:()=>({matches:false,addEventListener(){}}),getComputedStyle:()=>({getPropertyValue:()=> '#888'}),devicePixelRatio:1,addEventListener(){},setTimeout(){}};
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
