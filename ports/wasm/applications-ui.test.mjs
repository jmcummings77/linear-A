import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import vm from 'node:vm';
import {test} from 'node:test';

const source=readFileSync(new URL('../../applications/app.mjs',import.meta.url),'utf8');
const template=readFileSync(new URL('../../applications/template.html',import.meta.url),'utf8');

function playground({width=340,height=310}={}){
 const nodes=new Map(),timers=new Map(),workers=[],events=new Map();let nextTimer=0;
 class CanvasContext{
  constructor(){this.texts=[];this.strokes=[];this.path=[];this.boxes=[];}
  finite(...values){assert.ok(values.every(Number.isFinite),`Nonfinite drawing coordinates: ${values}`);}
  setTransform(...values){this.finite(...values);}
  clearRect(){this.texts=[];this.strokes=[];this.boxes=[];}
  beginPath(){this.path=[];}
  moveTo(x,y){this.finite(x,y);this.path.push([x,y]);}
  lineTo(x,y){this.finite(x,y);this.path.push([x,y]);}
  arc(...values){this.finite(...values);}
  stroke(){this.strokes.push({path:this.path.slice(),width:this.lineWidth,alpha:this.globalAlpha,color:this.strokeStyle});}
  fill(){}
  setLineDash(){}
  fillText(text,x,y){this.finite(x,y);this.texts.push({text,x,y,color:this.fillStyle,align:this.textAlign});}
  measureText(text){return {width:text.length*7};}
  fillRect(x,y,w,h){this.finite(x,y,w,h);this.boxes.push({x,y,w,h});}
  createImageData(w,h){return {data:new Uint8ClampedArray(w*h*4)};}
  putImageData(){}
 }
 class Element{
  constructor(){this.children=[];this.style={};this.textContent='';this.className='';this.hidden=false;this.clientWidth=width;this.clientHeight=height;this.ctx=new CanvasContext();this._value='';}
  set value(v){this._value=String(v);}
  get value(){return this._value;}
  append(...children){this.children.push(...children);}
  replaceChildren(...children){this.children=children;}
  getContext(){return this.ctx;}
  getBoundingClientRect(){return {left:0,top:0,width:this.clientWidth,height:this.clientHeight};}
  setPointerCapture(){}
 }
 const get=id=>{if(!nodes.has(id))nodes.set(id,new Element());return nodes.get(id);};
 get('bundle').textContent=JSON.stringify({sha256:'saved-test-bundle',worker_source:''});
 for(const [id,value] of Object.entries({'image-rank':6,'image-ink':1,'image-preset':'shapes','fit-degree':2}))get(id).value=value;
 class Worker{
  constructor(){workers.push(this);}
  postMessage(message){this.message=message;}
  terminate(){}
 }
 const context={document:{getElementById:get,createElement:()=>new Element(),documentElement:{}},Worker,Blob,URL:{createObjectURL:()=> 'blob:test',revokeObjectURL(){}},getComputedStyle:()=>({getPropertyValue:name=>({'--ink':'#111','--muted':'#555','--green':'#070','--orange':'#b50','--plot':'#eee'}[name])}),devicePixelRatio:2,setTimeout:(callback,delay)=>{timers.set(++nextTimer,{callback,delay});return nextTimer;},clearTimeout:id=>timers.delete(id),addEventListener:(name,callback)=>events.set(name,callback)};
 vm.runInNewContext(source,context);
 function flush(){for(const [id,{callback,delay}] of [...timers])if(delay===300){timers.delete(id);callback();}}
 function deliver(kind,result){flush();const worker=workers.find(w=>w.message?.kind===kind);assert.ok(worker,`Missing ${kind} worker`);worker.onmessage({data:{id:worker.message.id,result}});}
 return {get,deliver,flush,resize:()=>events.get('resize')()};
}
const pcaResult=({ambiguous=false,mean=[0,0],values=[25,.1],axes=[[1,0],[0,1]]}={})=>({mean,values,axes,explained:ambiguous?null:.9,ambiguous,projected:[]});
const imageResult=(rank,values=Array.from({length:24},(_,i)=>24-i))=>({singularValues:values,reconstructed:Array(576).fill(0),relativeError:0,rmse:0,retained:rank/24,factorScalars:rank*49,originalScalars:576});

test('both coordinate plots show sparse ticks and stronger zero axes after point selection',()=>{
 const page=playground();
 for(const kind of ['pca','fit']){
  const canvas=page.get(kind+'-canvas'),ctx=canvas.ctx;
  const labels=ctx.texts.map(item=>item.text);
  for(const expected of ['x','y','0','-4','-2','2','4'])assert.ok(labels.includes(expected),`${kind} misses ${expected}`);
  assert.equal(labels.filter(label=>label==='-4').length,2);assert.ok(!labels.includes('1'));
  const zero=ctx.strokes.find(stroke=>stroke.path.length===4&&stroke.path[0][1]===canvas.clientHeight/2&&stroke.path[2][0]===canvas.clientWidth/2);
  const grid=ctx.strokes.find(stroke=>stroke.alpha<1);assert.ok(zero&&grid);assert.ok(zero.width>grid.width);assert.ok(zero.alpha>grid.alpha);
  page.get(kind+'-selected').value=1;page.get(kind+'-selected').onchange();
  assert.equal(page.get(kind+'-x').value,String(-4+8/11));assert.ok(ctx.texts.some(label=>label.text==='x'));
 }
});

test('PCA labels remain within narrow canvas edges and disappear for ambiguous directions',()=>{
 const page=playground({width:280,height:310}),canvas=page.get('pca-canvas');
 page.deliver('pca',pcaResult({mean:[4.8,4.8]}));
 let labels=canvas.ctx.texts.filter(label=>/^PC[12]$/.test(label.text));assert.deepEqual(labels.map(label=>label.text),['PC1','PC2']);
 for(const label of labels){assert.ok(label.x>=0&&label.x+21<=canvas.clientWidth);assert.ok(label.y>=13&&label.y<=canvas.clientHeight);}
 assert.equal(canvas.ctx.boxes.length,2);
 const [a,b]=canvas.ctx.boxes;assert.ok(a.x+a.w<=b.x||b.x+b.w<=a.x||a.y+a.h<=b.y||b.y+b.h<=a.y,'PC tags overlap');
 page.get('pca-reset').onclick();page.deliver('pca',pcaResult({ambiguous:true,values:[0,0]}));
 assert.match(page.get('pca-status').textContent,/No unique principal direction/);
 assert.equal(canvas.ctx.texts.filter(label=>/^PC[12]$/.test(label.text)).length,0);
 assert.ok(canvas.ctx.texts.some(label=>label.text==='x'));
});

test('rank selection, cutoff and accessible values agree at zero, partial and full rank',()=>{
 const page=playground();
 for(const rank of [0,6,24]){
  page.get('image-rank').value=rank;page.get('image-rank').oninput();
  assert.equal(page.get('image-rank-label').textContent,String(rank));
  assert.equal(page.get('image-spectrum-cutoff').hidden,true);
  page.deliver('image',imageResult(rank));
  const bars=page.get('image-spectrum').children,rows=page.get('image-spectrum-values').children;
  assert.equal(bars.length,24);assert.equal(bars.filter(bar=>bar.className==='retained').length,rank);
  assert.equal(rows.filter(row=>row.children[2].textContent==='Retained').length,rank);
  assert.equal(page.get('image-spectrum-cutoff').hidden,false);
  assert.equal(page.get('image-spectrum-cutoff-label').textContent,`k = ${rank}`);
  if(rank===0){assert.equal(page.get('image-spectrum-cutoff').style.left,'0%');assert.match(page.get('image-spectrum-caption').textContent,/all 24 components omitted/);}
  else if(rank===24){assert.equal(page.get('image-spectrum-cutoff').style.left,'100%');assert.match(page.get('image-spectrum-caption').textContent,/all 24 components retained/);}
  else{assert.match(page.get('image-spectrum-caption').textContent,/retain σ1–σ6; omit σ7–σ24/);assert.match(page.get('image-spectrum-cutoff').style.left,/25%/);}
 }
 const rankControl=template.match(/<input\b[^>]*\bid="image-rank"[^>]*>/)?.[0]||'';
 assert.match(rankControl,/aria-label="Retained rank"/);
 assert.match(rankControl,/aria-describedby="image-spectrum-caption image-spectrum-legend"/);
});

test('blank-image spectra stay finite and editing clears outdated spectrum evidence',()=>{
 const page=playground();page.deliver('image',imageResult(6,Array(24).fill(0)));
 assert.ok(page.get('image-spectrum').children.every(bar=>bar.style.height==='0%'));
 page.get('image-preset').value='blank';page.get('image-preset').onchange();
 assert.equal(page.get('image-spectrum').children.length,0);assert.equal(page.get('image-spectrum-values').children.length,0);
 assert.equal(page.get('image-spectrum-cutoff').hidden,true);assert.match(page.get('image-spectrum-caption').textContent,/Calculating singular values/);
});
