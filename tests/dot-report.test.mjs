import assert from 'node:assert/strict';
import test from 'node:test';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';
const source=readFileSync(new URL('../experiments/machine-code-dot/report.mjs',import.meta.url),'utf8');
class Element {
 constructor(tag='div'){this.tagName=tag;this.children=[];this.style={};this.className='';this.value='';this._text='';}
 append(node){this.children.push(node);}
 replaceChildren(...nodes){this.children=nodes;this._text='';}
 set textContent(value){this._text=String(value);this.children=[];}
 get textContent(){return this._text+this.children.map(node=>node.textContent).join('');}
}
const timing=(implementation,median_ns,n=256,samples=[median_ns*.9,median_ns,median_ns*1.1])=>({implementation,median_ns,n,ns_per_call:samples});
function report(timings){
 const data={metadata:{timestamp_utc:'2026-10-02T12:00:00Z'},samples:3,correctness_checks:9,words:['0x1234'],timings};
 const original=JSON.stringify(data),nodes=new Map(),get=id=>{if(!nodes.has(id))nodes.set(id,new Element());return nodes.get(id);};
 get('data').textContent=original;get('names').textContent=JSON.stringify({'machine-code':'Hand-encoded machine code',fast:'Fast',slow:'Slow',slower:'Slower',zero:'Zero'});
 vm.runInNewContext(source,{document:{getElementById:get,createElement:tag=>new Element(tag)}});
 return {get,data,original,rows:()=>get('latency-chart').children.filter(row=>row.className.split(' ').includes('latency-row')&&!row.className.includes('latency-axis')),choose(n){get('size').value=n;get('size').onchange();}};
}
const child=(node,cls)=>node.children.find(n=>n.className.split(' ').includes(cls));
test('shared logarithmic positions preserve multiplicative distance and align the machine-code reference',()=>{
 const page=report([timing('slower',10000),timing('machine-code',100),timing('slow',1000),timing('fast',10)]);
 const rows=page.rows(),dots=rows.map(row=>parseFloat(child(child(row,'latency-track'),'latency-dot').style.left));
 assert.deepEqual(rows.map(row=>child(row,'latency-label').children[0].textContent),['Fast','Hand-encoded machine code','Slow','Slower']);
 assert.ok(Math.abs((dots[1]-dots[0])-(dots[2]-dots[1]))<1e-9);
 assert.ok(Math.abs((dots[2]-dots[1])-(dots[3]-dots[2]))<1e-9);
 for(const row of rows)assert.equal(parseFloat(child(child(row,'latency-track'),'latency-reference').style.left),dots[1]);
 const ticks=child(page.get('latency-chart').children[0],'latency-track').children.map(n=>n.textContent);
 assert.ok(ticks.includes('1×'));assert.ok(ticks.includes('10×'));assert.ok(ticks.includes('100×'));
 assert.match(page.get('latency-note').textContent,/Whiskers.*sample ranges/);
 assert.equal(JSON.stringify(page.data),page.original);
});
test('sample whiskers fit the domain, and changing vector length replaces the chart and exact table',()=>{
 const page=report([timing('machine-code',100),timing('slow',1000,256,[10,1000,100000]),timing('machine-code',5,16),timing('fast',2.5,16)]);
 page.choose(256);
 for(const row of page.rows()){
  const range=child(child(row,'latency-track'),'latency-range');
  assert.ok(parseFloat(range.style.left)>=0);
  assert.ok(parseFloat(range.style.left)+parseFloat(range.style.width)<=100+1e-9);
 }
 page.choose(16);
 assert.equal(page.rows().length,2);assert.equal(page.get('rows').children.length,2);
 assert.match(page.get('rows').children[0].textContent,/Fast.*2.5 ns.*0.50×/);
 assert.equal(JSON.stringify(page.data),page.original);
});
test('equal values and zero samples remain finite; zero medians are explicitly outside the logarithmic scale',()=>{
 const page=report([timing('machine-code',10,256,[10,10,10]),timing('fast',10,256,[10,10,10])]);
 for(const row of page.rows())assert.equal(child(child(row,'latency-track'),'latency-dot').style.left,'50%');
 const zero=report([timing('machine-code',10),timing('zero',0,256,[0,0,0])]);
 const zeroTrack=child(zero.rows()[0],'latency-track');
 assert.equal(child(zeroTrack,'latency-dot'),undefined);
 assert.match(zeroTrack.textContent,/outside log scale/);
 assert.match(zero.get('latency-note').textContent,/Zero timings are excluded/);
 assert.match(zero.get('rows').children[0].textContent,/0.0 ns/);
});
test('a zero or absent reference keeps exact measurements and avoids undefined ratios',()=>{
 for(const timings of [[timing('machine-code',0,256,[0,0,0]),timing('slow',10)],[timing('slow',10)]]){
  const page=report(timings);
  assert.equal(page.get('latency-chart').children.length,0);
  assert.match(page.get('latency-note').textContent,/positive machine-code timing/);
  assert.equal(page.get('rows').children.length,timings.length);
  assert.ok(page.get('rows').children.every(row=>row.children[3].textContent==='—'));
 }
});
